"""foilgeom: pure-Python geometry core for the FoilBuilder Fusion add-in.

No numpy (Fusion's bundled Python doesn't ship it). All lengths in mm.
Axes: x aft (chordwise), y spanwise, z up. Origin on the planform axis at the
centreline, on the span line.

Planform: bi-elliptical. LE and TE are semi-ellipses sharing a spanwise axis
placed DELTA * c_root ahead of the root TE:
    c(y) = c_root * sqrt(1 - (2y/b)^2),    c_root = 4 A / (pi b)
Anhedral (optional, "banana" foil): span line is an elliptical arc
    z(y) = -E_V * (1 - sqrt(1 - (y/E_H)^2)),  E_H solved for the tip drop.
Sections are rotated about x to stay perpendicular to the span line.
Flap (optional): straight hinge line at FLAP_FRAC * c_root ahead of the root
TE, ending where the local flap chord drops to FLAP_MIN_CHORD.
Sections (optional): SectionedAirfoil blends several airfoils placed at
fractions of the semi-span, so the section shape can change along the span.
"""
import math

# ------------------------------------------------------------------ airfoil I/O


def _floats(line):
    out = []
    for tok in line.replace(",", " ").split():
        try:
            out.append(float(tok))
        except ValueError:
            return None
    return out


def parse_dat(text):
    """Return (name, pts) with pts in Selig order (TE -> upper -> LE -> lower -> TE),
    normalised to unit chord. Accepts Selig and Lednicer formats."""
    lines = [l.strip() for l in text.splitlines()]
    name = "airfoil"
    rows = []
    for i, l in enumerate(lines):
        if not l:
            rows.append(None)
            continue
        f = _floats(l)
        if f is None or len(f) < 2:
            if i == 0 or not rows:
                name = l
            continue
        rows.append((f[0], f[1]))
    # Lednicer: first numeric row holds point counts (>1)
    num = [r for r in rows if r is not None]
    if not num:
        raise ValueError("No coordinate rows found in .dat file")
    if num[0][0] > 1.5 and num[0][1] > 1.5:
        n_up, n_lo = int(round(num[0][0])), int(round(num[0][1]))
        coords = num[1:]
        if len(coords) < n_up + n_lo:
            raise ValueError("Lednicer file: fewer points than declared")
        upper = coords[:n_up]            # LE -> TE
        lower = coords[n_up:n_up + n_lo]  # LE -> TE
        pts = list(reversed(upper)) + lower[1:]
    else:
        pts = num
    # drop consecutive duplicates
    clean = [pts[0]]
    for p in pts[1:]:
        if abs(p[0] - clean[-1][0]) > 1e-9 or abs(p[1] - clean[-1][1]) > 1e-9:
            clean.append(p)
    pts = clean
    xs = [p[0] for p in pts]
    xmin, xmax = min(xs), max(xs)
    c = xmax - xmin
    if c <= 0:
        raise ValueError("Degenerate airfoil (zero chord)")
    pts = [((x - xmin) / c, y / c) for x, y in pts]
    # orientation check: upper surface first (Selig). If the first half is
    # below the second half, reverse.
    i_le = min(range(len(pts)), key=lambda i: pts[i][0])
    up = pts[:i_le + 1]
    lo = pts[i_le:]
    if sum(p[1] for p in up) / len(up) < sum(p[1] for p in lo) / len(lo):
        pts = list(reversed(pts))
    # A closed TE (first point == last point) is fine: the TE opening applied in
    # Planform.section() separates the two points.
    return name, pts


def _interp(x, xp, fp):
    if x <= xp[0]:
        return fp[0]
    if x >= xp[-1]:
        return fp[-1]
    lo, hi = 0, len(xp) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xp[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - xp[lo]) / (xp[hi] - xp[lo]) if xp[hi] != xp[lo] else 0.0
    return fp[lo] + t * (fp[hi] - fp[lo])


def airfoil_stats(pts):
    """(t/c, x of max t, max camber, x of max camber) for unit-chord Selig points."""
    i_le = min(range(len(pts)), key=lambda i: pts[i][0])
    up = sorted(pts[:i_le + 1])
    lo = sorted(pts[i_le:])
    ux, uy = [p[0] for p in up], [p[1] for p in up]
    lx, ly = [p[0] for p in lo], [p[1] for p in lo]
    best_t = (0, 0)
    best_c = (0, 0)
    n = 400
    for k in range(n + 1):
        x = k / n
        yu, yl = _interp(x, ux, uy), _interp(x, lx, ly)
        t, cm = yu - yl, 0.5 * (yu + yl)
        if t > best_t[0]:
            best_t = (t, x)
        if cm > best_c[0]:
            best_c = (cm, x)
    return best_t[0], best_t[1], best_c[0], best_c[1]


# ------------------------------------------------------------------ spanwise sections

N_RESAMPLE = 50     # points per surface after resampling (cosine spaced)


def _cubic_spline(xs, ys):
    """Natural cubic spline through (xs, ys), xs strictly increasing. Returns f(x)."""
    n = len(xs)
    if n < 3:
        return lambda x: _interp(x, xs, ys)
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    a, b, c, d = [0.0] * n, [1.0] * n, [0.0] * n, [0.0] * n
    for i in range(1, n - 1):
        a[i], b[i], c[i] = h[i - 1], 2.0 * (h[i - 1] + h[i]), h[i]
        d[i] = 6.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])
    for i in range(1, n):                       # Thomas algorithm
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    m = [0.0] * n
    m[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        m[i] = (d[i] - c[i] * m[i + 1]) / b[i]

    def f(x):
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= x:
                lo = mid
            else:
                hi = mid
        x = min(max(x, xs[lo]), xs[hi])
        A = (xs[hi] - x) / h[lo]
        B = 1.0 - A
        return A * ys[lo] + B * ys[hi] + ((A ** 3 - A) * m[lo] + (B ** 3 - B) * m[hi]) * h[lo] ** 2 / 6.0
    return f


def resample(pts, n=N_RESAMPLE):
    """Resample unit-chord Selig points onto common cosine-spaced x stations, so
    that different airfoils can be blended point by point. Each surface is
    splined as y(sqrt(x)), which stays smooth round a rounded leading edge."""
    i_le = min(range(len(pts)), key=lambda i: pts[i][0])
    xs = [0.5 * (1.0 - math.cos(math.pi * k / n)) for k in range(n + 1)]
    surf = []
    for part in (pts[:i_le + 1], pts[i_le:]):
        us, ys = [], []
        for x, y in sorted(part):
            u = math.sqrt(max(x, 0.0))
            if us and u - us[-1] < 1e-9:
                continue
            us.append(u)
            ys.append(y)
        f = _cubic_spline(us, ys)
        surf.append([(x, f(math.sqrt(x))) for x in xs])
    up, lo = surf
    return list(reversed(up)) + lo[1:]


class SectionedAirfoil:
    """Airfoil that changes along the span.

    sections: [(eta, name, pts)], eta = |y| / semi-span (0 centreline, 1 tip),
    mirrored to both sides. Between two stations the resampled coordinates are
    blended linearly; inboard of the first / outboard of the last station that
    section is used as it is."""

    def __init__(self, sections, n=N_RESAMPLE):
        if not sections:
            raise ValueError("Add at least one section")
        secs = sorted(sections, key=lambda s: s[0])
        for i, (eta, name, _) in enumerate(secs):
            if not 0.0 <= eta <= 1.0:
                raise ValueError(f"Section position {eta} of {name} must be between 0 and 1")
            if i and eta - secs[i - 1][0] < 1e-6:
                raise ValueError(f"Two sections at the same position ({eta:.3f})")
        self.etas = [s[0] for s in secs]
        self.names = [s[1] for s in secs]
        self.afs = [resample(s[2], n) for s in secs]
        self._cache = {}    # keeps blended lists alive (Planform caches stats by id)

    @property
    def name(self):
        out = []
        for nm in self.names:
            if nm not in out:
                out.append(nm)
        return " / ".join(out)

    def at(self, eta):
        eta = min(abs(eta), 1.0)
        key = round(eta, 9)
        if key in self._cache:
            return self._cache[key]
        e = self.etas
        if eta <= e[0]:
            af = self.afs[0]
        elif eta >= e[-1]:
            af = self.afs[-1]
        else:
            i = max(k for k in range(len(e)) if e[k] <= eta)
            t = (eta - e[i]) / (e[i + 1] - e[i])
            af = [(xa, (1.0 - t) * ya + t * yb)
                  for (xa, ya), (_, yb) in zip(self.afs[i], self.afs[i + 1])]
        self._cache[key] = af
        return af


# ------------------------------------------------------------------ parameters

TE_MIN_MM = 0.05

PRESETS = {
    "Bonetti main foil (E393)": dict(span=1150.0, area_cm2=960.0, delta=0.43, flap_on=True,
                                     flap_frac=0.30, flap_min_chord=10.0, anhedral_depth=0.0,
                                     anhedral_semiaxis=65.0, n_half=18, tip_eta=0.998,
                                     te_thick_pct=0.4, bulb_on=True, bulb_length=200.0,
                                     bulb_nose_frac=0.6, bulb_thickness=24.0, bulb_half_width=30.0),
    "Bonetti rudder foil (E393)": dict(span=800.0, area_cm2=400.0, delta=0.30, flap_on=False,
                                       flap_frac=0.30, flap_min_chord=10.0, anhedral_depth=40.0,
                                       anhedral_semiaxis=65.0, n_half=18, tip_eta=0.998,
                                       te_thick_pct=0.4, bulb_on=True, bulb_length=170.0,
                                       bulb_nose_frac=0.6, bulb_thickness=17.0, bulb_half_width=25.0),
}
STRUT_DEFAULTS = dict(strut_on=False, strut_chord=80.0, strut_le_x=15.0, strut_height=150.0,
                      strut_rim_round=1.0, strut_slot_depth=5.0)
STRUT_MIN_FLOOR_MM = 2.0      # foil material left under the slot
STRUT_TE_MIN_MM = 0.5         # vertical's TE thickness (rounded off in CAD)
DEFAULTS = dict(PRESETS["Bonetti main foil (E393)"])
DEFAULTS.update(STRUT_DEFAULTS)


class Planform:
    def __init__(self, p):
        self.p = p
        self.span = float(p["span"])
        self.area = float(p["area_cm2"]) * 100.0          # mm^2
        self.delta = float(p["delta"])
        self.b2 = self.span / 2.0
        self.c_root = 4.0 * self.area / (math.pi * self.span)
        d, ev = float(p["anhedral_depth"]), float(p["anhedral_semiaxis"])
        if d > 0:
            if d >= ev:
                raise ValueError("Anhedral depth must be smaller than the ellipse semi-axis")
            r = 1.0 - d / ev
            self.e_h = self.b2 / math.sqrt(1.0 - r * r)
            self.e_v = ev
        else:
            self.e_h = None
            self.e_v = 0.0
        # blended bulb (optional; older parameter sets without bulb keys -> no bulb)
        self.bulb_on = bool(p.get("bulb_on", False))
        self.bulb_w = float(p.get("bulb_half_width", 30.0))
        self.bulb_nose = float(p.get("bulb_nose_frac", 0.6))
        self.bulb_t = float(p.get("bulb_thickness", 0.0))
        self.bulb_dL = max(0.0, float(p.get("bulb_length", 0.0)) - self.c_root) if self.bulb_on else 0.0
        if self.bulb_on:
            if not (0 < self.bulb_w < 0.4 * self.b2):
                raise ValueError("Bulb half-width must be between 0 and 40% of the semi-span")
            if not (0.0 <= self.bulb_nose <= 1.0):
                raise ValueError("Bulb nose fraction must be between 0 and 1")
            if self.bulb_t <= 0:
                raise ValueError("Bulb thickness must be positive")
        self._taf = {}

    # span line
    def z_span(self, y):
        if self.e_h is None:
            return 0.0, 0.0
        q = math.sqrt(max(1.0 - (y / self.e_h) ** 2, 1e-12))
        return -self.e_v * (1.0 - q), -self.e_v * (y / self.e_h ** 2) / q

    def chord(self, y):
        return self.c_root * math.sqrt(max(1.0 - (y / self.b2) ** 2, 0.0))

    def stations(self, extra_eta=()):
        """Span stations (mm). extra_eta: semi-span fractions that must get a
        station of their own (e.g. where a user-defined section sits); they
        replace any regular station closer than 1% of the semi-span."""
        n = int(self.p["n_half"])
        eta_tip = float(self.p["tip_eta"])
        out = []
        for k in range(-n, n + 1):
            th = (math.pi / 2.0) * k / n
            out.append(math.sin(th) * eta_tip * self.b2)
        y_in = 0.0
        if self.bulb_on:
            yb = 1.15 * self.bulb_w
            nb = 8
            y_in = yb * 1.05
            out = [y for y in out if abs(y) > y_in]
            out += [yb * k / nb for k in range(-nb, nb + 1)]
        gap = 0.01 * self.b2
        for e in extra_eta:
            ye = abs(e) * self.b2
            if ye <= max(y_in, gap) or ye >= eta_tip * self.b2 - gap:
                continue            # already covered by centreline / bulb / tip stations
            out = [y for y in out if abs(abs(y) - ye) > gap]
            out += [-ye, ye]
        out.sort()
        return out

    # ---- blended bulb
    def bulb_f(self, y):
        """Blend weight: 1 on the centreline, 0 beyond bulb_half_width (C1 smooth)."""
        if not self.bulb_on or abs(y) >= self.bulb_w:
            return 0.0
        return math.cos(0.5 * math.pi * y / self.bulb_w) ** 2

    def bulb_tail_x(self):
        return self.delta * self.c_root + self.bulb_dL * (1.0 - self.bulb_nose)

    def _t_af(self, af):
        k = id(af)
        if k not in self._taf:
            self._taf[k] = airfoil_stats(af)[0]
        return self._taf[k]

    def flap(self):
        """(x_hinge, y_flap_half) or None."""
        if not self.p.get("flap_on"):
            return None
        x_te_root = self.delta * self.c_root
        x_h = x_te_root - float(self.p["flap_frac"]) * self.c_root
        s_min = (x_h + float(self.p["flap_min_chord"])) / x_te_root
        if s_min >= 1.0:
            return None
        return x_h, self.b2 * math.sqrt(1.0 - s_min * s_min)

    def flap_cutter_polygon(self):
        """Planform polygon (x, y) in mm of the flap cutter, or None.
        With a bulb, the cutter goes round the bulb tail, giving port and
        starboard flap halves that start where the bulb has blended out."""
        f = self.flap()
        if not f:
            return None
        x_h, y_f = f
        x_far = self.delta * self.c_root + self.bulb_dL + 300.0
        if not self.bulb_on:
            return [(x_h, -y_f), (x_far, -y_f), (x_far, y_f), (x_h, y_f)]
        y_in = self.bulb_w
        if y_in >= y_f:
            return None
        x_t = self.bulb_tail_x() + 5.0
        return [(x_h, -y_f), (x_far, -y_f), (x_far, y_f), (x_h, y_f),
                (x_h, y_in), (x_t, y_in), (x_t, -y_in), (x_h, -y_in)]

    def section(self, af, y):
        """3D points (mm) of the section at span station y, Selig order."""
        c_w = self.chord(y)
        fb = self.bulb_f(y)
        d_fwd = self.bulb_dL * self.bulb_nose * fb
        d_aft = self.bulb_dL * (1.0 - self.bulb_nose) * fb
        c = c_w + d_fwd + d_aft                      # local chord incl. bulb
        zs = c                                       # vertical scale
        if fb > 0:
            k0 = self.bulb_t / (self._t_af(af) * (self.c_root + self.bulb_dL))
            zs = c * (1.0 + (k0 - 1.0) * fb)
        # TE thickness as a fraction of chord, never below TE_MIN_MM (avoids
        # microscopic edges at the tip that CAD kernels dislike)
        te = max(float(self.p["te_thick_pct"]) / 100.0, TE_MIN_MM / max(c, 1e-9))
        i_le = min(range(len(af)), key=lambda i: af[i][0])
        x_le = -(1.0 - self.delta) * c_w - d_fwd
        z0, dz = self.z_span(y)
        phi = math.atan(dz)
        ny, nz = -math.sin(phi), math.cos(phi)
        pts = []
        for i, (xa, za) in enumerate(af):
            sgn = 1.0 if i <= i_le else -1.0
            if i == i_le:
                sgn = 0.0
            zz = za * zs + sgn * xa * te / 2.0 * c
            pts.append((x_le + xa * c, y + zz * ny, z0 + zz * nz))
        return pts

    def summary(self):
        f = self.flap()
        s = {"c_root": self.c_root, "AR": self.span ** 2 / self.area,
             "x_le_root": -(1 - self.delta) * self.c_root, "x_te_root": self.delta * self.c_root}
        if self.e_h is not None:
            y_t = float(self.p["tip_eta"]) * self.b2
            s["tip_drop"] = -self.z_span(self.b2)[0]
            s["tip_slope_deg"] = math.degrees(math.atan(abs(self.z_span(y_t)[1])))
        if f:
            s["x_hinge"], s["flap_span"] = f[0], 2 * f[1]
        if self.bulb_on:
            s["bulb_length"] = self.c_root + self.bulb_dL
            s["bulb_nose_ahead_of_le"] = self.bulb_dL * self.bulb_nose
            s["bulb_tail_x"] = self.bulb_tail_x()
            s["bulb_thickness"] = self.bulb_t
        return s


# ------------------------------------------------------------------ vertical support


def _surfaces(pts):
    """3D section points (x, y, z), Selig order -> upper and lower surface as
    x-sorted (xs, zs) lists."""
    i_le = min(range(len(pts)), key=lambda i: pts[i][0])
    out = []
    for part in (pts[:i_le + 1], pts[i_le:]):
        s = sorted((q[0], q[2]) for q in part)
        out.append(([q[0] for q in s], [q[1] for q in s]))
    return out


class Strut:
    """Vertical support leaving the top of the foil at the centreline.

    A straight, separate piece. Foil axes: chord along x, thickness along y,
    straight up along z to strut_height (measured from the foil chord line at
    the centreline), where it ends in a flat top face. At the bottom it reaches
    strut_slot_depth below the lowest point of the foil top under it; the same
    shape is cut out of the foil as the slot, whose rim can be rounded by
    strut_rim_round (done in CAD).
    eta = 0 at the highest point of the foil top under the vertical, 1 at its top."""

    def __init__(self, p, pf, foil_af_for, strut_afs, n=200):
        self.c = float(p["strut_chord"])
        self.h = float(p["strut_height"])
        self.rim = float(p.get("strut_rim_round", 0.0))
        self.slot = float(p["strut_slot_depth"])
        if self.c <= 0 or self.slot <= 0 or self.rim < 0:
            raise ValueError("Vertical chord and slot depth must be positive, the rim round 0 or more")
        # strut_le_x: LE position aft of the foil nose (bulb nose with a bulb) on the centreline
        self.x_le = _surfaces(pf.section(foil_af_for(0.0), 0.0))[0][0][0] + float(p["strut_le_x"])
        self.x_te = self.x_le + self.c
        # the CAD side closes this section with a smooth rounded TE, which needs a
        # little room to turn: never thinner than STRUT_TE_MIN_MM
        self.te = max(float(p["te_thick_pct"]) / 100.0, STRUT_TE_MIN_MM / self.c)
        self.half_t = max(airfoil_stats(a)[0] for a in strut_afs) * self.c / 2.0
        margin = self.rim + 1.0               # foil material kept round the slot

        # foil surfaces under the vertical, on the centreline and at its sides
        top_min, low_max = float("inf"), -float("inf")
        for y in (0.0, self.half_t):
            (ux, uz), (lx, lz) = _surfaces(pf.section(foil_af_for(y), y))
            if y == 0.0 and (self.x_le - margin <= ux[0] or self.x_te + margin >= ux[-1]):
                raise ValueError("The vertical must sit within the foil chord at the centreline, "
                                 f"with {margin:g} mm to spare (check chord / LE position)")
            xs = [self.x_le + self.c * k / n for k in range(n + 1)]
            if y == 0.0:
                self.z_root = max(_interp(x, ux, uz) for x in xs)
            top_min = min(top_min, min(_interp(x, ux, uz) for x in xs))
            low_max = max(low_max, max(_interp(x, lx, lz) for x in xs))
        self.z_bot = top_min - self.slot
        self.floor = self.z_bot - low_max
        if self.floor < STRUT_MIN_FLOOR_MM:
            raise ValueError(f"Slot too deep: only {self.floor:.1f} mm of foil left under it "
                             f"(min {STRUT_MIN_FLOOR_MM:.0f} mm)")
        h_min = self.z_root + 5.0
        if self.h < h_min:
            raise ValueError(f"Vertical height must be at least {h_min:.0f} mm (foil top + 5 mm)")
        if pf.flap_cutter_polygon():
            if pf.bulb_on:
                bad = (self.half_t + margin >= pf.bulb_w
                       or self.x_te + margin >= pf.bulb_tail_x() + 5.0)
            else:
                bad = self.x_te + margin >= pf.flap()[0]
            if bad:
                raise ValueError("The vertical's slot overlaps the flap: move it forward, "
                                 "shorten its chord or turn the flap off")

    def eta(self, z):
        return min(max((z - self.z_root) / (self.h - self.z_root), 0.0), 1.0)

    def stations(self, etas=()):
        """Section heights (mm): tenon bottom and top; with sections also the
        foil top, each section's eta and 3 levels in between."""
        zs = {self.z_bot, self.h}
        if etas:
            pts = sorted(set([0.0, 1.0] + [min(max(e, 0.0), 1.0) for e in etas]))
            for a, b in zip(pts, pts[1:]):
                for k in range(4):
                    zs.add(self.z_root + (a + (b - a) * k / 4.0) * (self.h - self.z_root))
        return sorted(zs)

    def section(self, af, z):
        """3D points (mm) of the horizontal section at height z, Selig order."""
        i_le = min(range(len(af)), key=lambda i: af[i][0])
        pts = []
        for i, (xa, ya) in enumerate(af):
            sgn = 0.0 if i == i_le else (1.0 if i < i_le else -1.0)
            pts.append((self.x_le + xa * self.c, (ya + sgn * xa * self.te / 2.0) * self.c, z))
        return pts

    def summary(self):
        return {"strut_t": 2 * self.half_t, "strut_floor": self.floor,
                "strut_z_root": self.z_root, "strut_z_bot": self.z_bot}
