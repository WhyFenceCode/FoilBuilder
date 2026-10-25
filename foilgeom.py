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
STRUT_DEFAULTS = dict(strut_on=False, strut_chord=80.0, strut_le_x=25.0, strut_height=150.0,
                      strut_rim_round=0.5, strut_slot_depth=5.0)
STRUT_MIN_FLOOR_MM = 2.0      # foil material left under the slot
STRUT_TE_MIN_MM = 0.5         # vertical's TE thickness (rounded off in CAD)
BLEND_DEFAULTS = dict(blend_on=False, blend_height=3.0, blend_width=6.0, blend_front_width=6.0,
                      blend_back_width=6.0, blend_wall_angle=25.0)
BLEND_END_POWER = 6           # front/back lengths apply along +-x, falling off as cos^this sideways
BLEND_FOOT_DEG = 8.0          # the blend meets the foil this far off tangent (clean join)
BLEND_SKIRT_MM = 2.0          # ...then dips under the foil top for this far
BLEND_DEPTH_MM = 1.2          # the blend body reaches this far under the foil top (overlap to join)
BLEND_SINK_MM = 0.3           # where faded out, its surface sits this far under the foil top
BLEND_FADE_MM = (3.0, 9.0)    # it fades in as the foil thickens from / to (mm)
BLEND_MIN_THICK_MM = 1.5      # no blend body where the foil is thinner than this
BLEND_SLICE_PTS = (40, 12)    # points on each slice's top / bottom, evenly spaced along it
BLEND_SIDE_STEP_MM = 4.0      # slices this far apart along the vertical's sides...
BLEND_FAN_STEP_DEG = 12.0     # ...and every this many degrees round its nose and tail
BLEND_INSIDE_MM = 0.15        # slices start this far inside the wall (cut away by the slot)
BLEND_END_FADE_MM = 15.0      # an open ring (overhanging vertical) fades out over this
DEFAULTS = dict(PRESETS["Bonetti main foil (E393)"])
DEFAULTS.update(STRUT_DEFAULTS)
DEFAULTS.update(BLEND_DEFAULTS)


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
    the centreline), where it ends in a flat top face. It can sit anywhere
    along the foil, also overhanging its nose or tail.

    Its bottom is the tenon line (tenon(), in the x-z plane): over the foil it
    runs strut_slot_depth below the foil top, but never closer than
    STRUT_MIN_FLOOR_MM to the foil bottom (where the foil is thinner it rises
    to the surface); past the foil's nose or tail it is the height of that
    end, so the overhang is trimmed level with it instead of hanging below.
    The CAD side trims the vertical to this line and cuts it out of the foil
    as the slot, whose rim can be rounded by strut_rim_round.
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

        # foil top / bottom on the centreline and across the vertical's thickness
        self._rows = []
        for y in (0.0, 0.5 * self.half_t, self.half_t):
            (ux, uz), (lx, lz) = _surfaces(pf.section(foil_af_for(y), y))
            self._rows.append((_row_spline(ux, uz), _row_spline(lx, lz)))
        (top0, le0, te0), (low0, _, _) = self._rows[0]
        if self.x_te <= le0 + 1.0 or self.x_le >= te0 - 1.0:
            raise ValueError("The vertical must sit at least partly over the foil (check chord / LE position)")
        xs = [self.x_le + self.c * k / n for k in range(n + 1)]
        on = [x for x in xs if le0 <= x <= te0]
        self.z_root = max(top0(x) for x in on)
        self.overhang = (max(0.0, le0 - self.x_le), max(0.0, self.x_te - te0))

        # tenon line, sampled ~every mm a little beyond both ends of the vertical
        m = max(int(self.c + 4.0), 8)
        xs_t = [self.x_le - 2.0 + (self.c + 4.0) * k / m for k in range(m + 1)]
        self.tenon_pts = [(x, self._tenon_z(x)) for x in xs_t]
        self.z_bot = min(z for _, z in self.tenon_pts) - 1.0       # loft bottom, below the line
        # deepest slot (on the centreline), and the least foil left under it on
        # any slice where there is a slot at all. Normally >= STRUT_MIN_FLOOR_MM;
        # less only right at a blunt nose / tail the vertical overhangs, where
        # the line has to come down to that end's height to leave no gap.
        self.slot_max = max(top0(x) - self._tenon_z(x) for x in on)
        self.floor = min((z + 0.1 - low(x)
                          for x, z in self.tenon_pts
                          for (top, le, te), (low, _, _) in self._rows
                          if le <= x <= te and z < top(x) - 0.2), default=0.0)

        h_min = self.z_root + 5.0
        if self.h < h_min:
            raise ValueError(f"Vertical height must be at least {h_min:.0f} mm (foil top + 5 mm)")
        if pf.flap_cutter_polygon():
            margin = self.rim + 1.0
            if pf.bulb_on:
                bad = (self.half_t + margin >= pf.bulb_w
                       or self.x_te + margin >= pf.bulb_tail_x() + 5.0)
            else:
                bad = self.x_te + margin >= pf.flap()[0]
            if bad:
                raise ValueError("The vertical's slot overlaps the flap: move it forward, "
                                 "shorten its chord or turn the flap off")

    def _surf_at(self, x):
        """(lowest foil top, highest foil bottom) at x across the vertical's
        thickness. Where a row is off the foil, its nose / tail height counts
        as both, so the line there lies level with that end."""
        tops, lows = [], []
        for (top, le, te), (low, _, _) in self._rows:
            if le <= x <= te:
                tops.append(top(x))
                lows.append(low(x))
            else:
                z = top(le) if x < le else top(te)
                tops.append(z)
                lows.append(z)
        return min(tops), max(lows)

    def _tenon_z(self, x):
        """Height of the vertical's bottom at x: slot_depth below the foil top
        but at least STRUT_MIN_FLOOR_MM above its bottom, and never above the
        top (no gap at the vertical's sides). 0.1 mm into the foil."""
        t, b = self._surf_at(x)
        return min(t, max(t - self.slot, b + STRUT_MIN_FLOOR_MM)) - 0.1

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
        """3D points (mm) of the horizontal section at height z: Selig order,
        then a few points round a semicircular TE, as one closed loop (the CAD
        side closes it with a closed spline)."""
        i_le = min(range(len(af)), key=lambda i: af[i][0])
        pts = []
        for i, (xa, ya) in enumerate(af):
            sgn = 0.0 if i == i_le else (1.0 if i < i_le else -1.0)
            pts.append((self.x_le + xa * self.c, (ya + sgn * xa * self.te / 2.0) * self.c, z))
        (x0, y0, _), (x1, y1, _) = pts[0], pts[-1]
        xm, ym, r = 0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * math.hypot(x0 - x1, y0 - y1)
        for th in (-60.0, -30.0, 0.0, 30.0, 60.0):          # lower TE -> aft -> upper TE
            t = math.radians(th)
            pts.append((xm + r * math.cos(t), ym + r * math.sin(t), z))
        return pts

    def summary(self):
        return {"strut_t": 2 * self.half_t, "strut_floor": self.floor, "strut_slot_max": self.slot_max,
                "strut_z_root": self.z_root, "strut_overhang_fwd": self.overhang[0],
                "strut_overhang_aft": self.overhang[1]}


def _row_spline(xs, zs):
    """Smooth z(x) along one surface of a section: (f, x_le, x_te). Splined in
    sqrt(x - x_le), which stays smooth round a rounded leading edge."""
    x0 = xs[0]
    us, vs = [], []
    for x, z in zip(xs, zs):
        u = math.sqrt(max(x - x0, 0.0))
        if us and u - us[-1] < 1e-9:
            continue
        us.append(u)
        vs.append(z)
    f = _cubic_spline(us, vs)
    return (lambda x: f(math.sqrt(max(x - x0, 0.0)))), x0, xs[-1]


def _smoothstep(t):
    t = min(max(t, 0.0), 1.0)
    return t * t * (3.0 - 2.0 * t)


class Blend:
    """Collar on the foil round the vertical: the foil top rises to meet the
    vertical closer to its angle before the slot is cut.

    Built like a fillet: slices at right angles to the vertical's wall, all
    round it (fanning out round its nose and tail), each going from just
    inside the wall outwards. Along a slice, at distance s from the wall, the
    collar stands blend_height above the foil top at the wall and falls to 0
    at the slice's width (by its direction: blend_front_width straight ahead,
    blend_back_width straight behind, blend_width sideways, between them by
    cos^BLEND_END_POWER of the angle), along a cubic that leaves the vertical
    at blend_wall_angle from vertical and meets the foil BLEND_FOOT_DEG off
    tangent; beyond that it dips under the foil top (so the join in CAD is a
    clean crossing, not two surfaces touching). Just inside the wall it keeps
    rising, levelling off, so the slot cut crosses it cleanly. It is scaled
    by a fade (0..1) that goes to 0 where the foil gets thin and near the
    flap; where faded out its surface sits BLEND_SINK_MM under the foil top,
    i.e. nothing shows. Where the wall itself is off the foil (an overhanging
    vertical) there are no slices, and the ring is left open with its ends
    faded out. So it fits wherever the vertical is.

    Every slice has the same make-up (wall, collar, foot, hidden part) and the
    same number of evenly spaced points, so a loft through them pairs up like
    with like and the surface stays smooth.

    CAD: each slice is a lens between that surface and BLEND_DEPTH_MM under
    the foil top; lofted round (closed when it goes all the way round) and
    joined to the foil."""

    def __init__(self, p, pf, foil_af_for, st, strut_af_at):
        self.st = st
        self.w = float(p["blend_width"])
        self.w_front = float(p.get("blend_front_width", self.w))
        self.w_back = float(p.get("blend_back_width", self.w))
        self.h = float(p["blend_height"])
        wall = float(p["blend_wall_angle"])
        if min(self.w, self.w_front, self.w_back) <= 0 or self.h <= 0 or not 0.0 <= wall < 80.0:
            raise ValueError("Blend widths and height must be positive, the wall angle 0-80 deg")
        if st.h < st.z_root + self.h + 5.0:
            raise ValueError(f"Vertical height must be at least {st.z_root + self.h + 5.0:.0f} mm "
                             "(blend + 5 mm)")
        self.tan_foot = math.tan(math.radians(BLEND_FOOT_DEG))
        self.wall_deg = wall
        self.cot_wall = 1.0 / math.tan(math.radians(max(wall, 1.0)))

        # the vertical's outline as a closed loop, with cumulative arc length
        P = [(x, y) for x, y, _ in st.section(strut_af_at(0.0), 0.0)]
        area = sum(P[i - 1][0] * P[i][1] - P[i][0] * P[i - 1][1] for i in range(len(P)))
        self._ccw = 1.0 if area > 0 else -1.0       # (ty, -tx) * this points out
        self._poly = P
        self._L = [0.0]
        for i in range(1, len(P) + 1):
            (x0, y0), (x1, y1) = P[i - 1], P[i % len(P)]
            self._L.append(self._L[-1] + math.hypot(x1 - x0, y1 - y0))

        # foil top and bottom at |y| stations round the vertical, each a smooth
        # spline in x (the blend meets the foil nearly tangent, so kinks in a
        # piecewise-linear surface would show as wobbles in its edge)
        side = max((self.w + c ** BLEND_END_POWER * (max(self.w_front, self.w_back) - self.w))
                   * math.sqrt(1.0 - c * c) for c in [k / 50.0 for k in range(51)])
        self.y_max = st.half_t + max(side, self.w) + BLEND_SKIRT_MM + 1.0
        self._dy = (self.y_max + 1.0) / 24.0
        self._top, self._low = [], []
        for k in range(-1, 27):
            (ux, uz), (lx, lz) = _surfaces(pf.section(foil_af_for(k * self._dy), k * self._dy))
            self._top.append(_row_spline(ux, uz))
            self._low.append(_row_spline(lx, lz))

        # where the flap is split off (the blend fades out before it)
        self._flap = None
        if pf.flap_cutter_polygon():
            x_h = pf.flap()[0]
            self._flap = (x_h, pf.bulb_w, pf.bulb_tail_x() + 5.0) if pf.bulb_on else (x_h, None, None)

        self._stations = self._place_stations()

    # ---- foil surfaces near the centreline (rows at y = -dy .. 26 dy; smooth in y
    #      by Catmull-Rom across rows)
    def _across(self, vals_at, y):
        t = min(abs(y) / self._dy, 24.999999) + 1.0
        k = int(t)
        f = t - k
        p0, p1, p2, p3 = (vals_at(i) for i in (k - 1, k, k + 1, k + 2))
        return 0.5 * (2 * p1 + (p2 - p0) * f + (2 * p0 - 5 * p1 + 4 * p2 - p3) * f * f
                      + (3 * p1 - p0 - 3 * p2 + p3) * f ** 3)

    def top(self, x, y):
        return self._across(lambda i: self._top[i][0](x), y)

    def low(self, x, y):
        return self._across(lambda i: self._low[i][0](x), y)

    def _on_foil(self, x, y):
        """Foil thickness at (x, y), 0 off the foil."""
        le = self._across(lambda i: self._top[i][1], y)
        te = self._across(lambda i: self._top[i][2], y)
        return self.top(x, y) - self.low(x, y) if le < x < te else 0.0

    # ---- the vertical's outline
    def _point(self, s):
        """Point on the outline at arc length s (wraps round)."""
        L, P = self._L, self._poly
        s %= L[-1]
        lo, hi = 0, len(L) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (mid, hi) if L[mid] <= s else (lo, mid)
        t = (s - L[lo]) / max(L[hi] - L[lo], 1e-12)
        (x0, y0), (x1, y1) = P[lo], P[hi % len(P)]
        return x0 + t * (x1 - x0), y0 + t * (y1 - y0)

    def _normal(self, s, ds=0.1):
        """Outward unit normal of the outline at arc length s (from a short
        chord, so the closely spaced points round its TE don't make it jump)."""
        (x0, y0), (x1, y1) = self._point(s - ds), self._point(s + ds)
        tx, ty = x1 - x0, y1 - y0
        L = math.hypot(tx, ty) or 1.0
        return self._ccw * ty / L, -self._ccw * tx / L

    def _place_stations(self):
        """Arc lengths of the slices round the outline: every BLEND_SIDE_STEP_MM
        along the sides, every BLEND_FAN_STEP_DEG of turning round the ends."""
        total = self._L[-1]
        n = max(int(total / 0.05), 200)
        meas, prev = [0.0], None
        for k in range(n + 1):
            nx, ny = self._normal(total * k / n)
            ang = math.atan2(ny, nx)
            if prev is not None:
                da = abs((ang - prev + math.pi) % (2 * math.pi) - math.pi)
                meas.append(meas[-1] + (total / n) / BLEND_SIDE_STEP_MM + da / math.radians(BLEND_FAN_STEP_DEG))
            prev = ang
        m = max(int(round(meas[-1])), 8)
        out, j = [], 0
        for i in range(m):
            target = meas[-1] * i / m
            while meas[j + 1] < target:
                j += 1
            t = (target - meas[j]) / max(meas[j + 1] - meas[j], 1e-12)
            out.append(total * (j + t) / n)
        return out

    # ---- shape
    def _width(self, nx):
        """Blend width where the outline's outward normal has x-component nx."""
        f = max(-nx, 0.0) ** BLEND_END_POWER
        b = max(nx, 0.0) ** BLEND_END_POWER
        return self.w + f * (self.w_front - self.w) + b * (self.w_back - self.w)

    def _bump(self, d, w):
        """Height above the foil top at signed distance d from the wall, for
        width w, before fading."""
        m0 = max(-w * self.cot_wall / self.h, -3.0)             # Hermite end slopes (in d / w)
        m1 = -self.tan_foot * w / self.h
        if d <= 0.0:
            # inside the vertical (cut away by the slot): keep rising, levelling
            # off smoothly, so there is no corner whose ripples reach the outside
            return self.h * (1.0 - m0 / w * (1.0 - math.exp(d)))
        if d >= w:
            return -self.tan_foot * (d - w)
        t = d / w
        return self.h * ((2 * t ** 3 - 3 * t ** 2 + 1) + (t ** 3 - 2 * t ** 2 + t) * m0
                         + (t ** 3 - t ** 2) * m1)

    def _fade(self, x, y, thick):
        lo, hi = BLEND_FADE_MM
        f = _smoothstep((thick - lo) / (hi - lo))
        if self._flap:
            x_h, y_b, x_t = self._flap
            room = x_h - x if y_b is None else min(max(x_h - x, y_b - abs(y)), x_t - x)
            f *= _smoothstep((room - 1.0) / 4.0)
        return f

    @staticmethod
    def _even(curve, n):
        """n points evenly spaced along a dense polyline of (s, z): a loft pairs
        up points by their fraction along each slice, so every slice gets the
        same number, spaced the same way."""
        L = [0.0]
        for (s0, z0), (s1, z1) in zip(curve, curve[1:]):
            L.append(L[-1] + math.hypot(s1 - s0, z1 - z0))
        out, j = [], 1
        for k in range(n):
            t_s = L[-1] * k / (n - 1)
            while j < len(L) - 1 and L[j] < t_s:
                j += 1
            t = (t_s - L[j - 1]) / max(L[j] - L[j - 1], 1e-12)
            (s0, z0), (s1, z1) = curve[j - 1], curve[j]
            out.append((s0 + t * (s1 - s0), z0 + t * (z1 - z0)))
        return out

    def _reach(self, x, y, nx, ny, w):
        """How far out from the wall point (x, y) the slice goes: its width +
        skirt, or less where the foil gets too thin. None if the wall point
        itself isn't on enough foil."""
        ok = lambda s: self._on_foil(x + nx * s, y + ny * s) >= BLEND_MIN_THICK_MM
        if not ok(0.0):
            return None
        s_max = w + BLEND_SKIRT_MM
        s = 0.0
        while s < s_max:
            s_next = min(s + 0.25, s_max)
            if not ok(s_next):
                a, b = s, s_next
                for _ in range(30):
                    m = 0.5 * (a + b)
                    a, b = (m, b) if ok(m) else (a, m)
                return a
            s = s_next
        return s_max

    def runs(self):
        """[(closed, [slice, ...])]: the slices to loft, in order round the
        vertical. One closed run when it goes all the way round, else open runs
        (their ends faded out). Each slice is (x, y, nx, ny, upper, lower):
        wall point, outward normal, upper points (wall outwards) and lower
        points (back to the wall), 3D (mm)."""
        sts = self._stations
        info = []
        for s in sts:
            (x, y), (nx, ny) = self._point(s), self._normal(s)
            w = self._width(nx)
            r = self._reach(x, y, nx, ny, w)
            info.append((s, x, y, nx, ny, w, r if r is not None and r >= 1.0 else None))
        valid = [i[6] is not None for i in info]
        if not any(valid):
            return []
        if all(valid):
            groups, closed = [list(range(len(info)))], True
        else:
            start = valid.index(False)                # walk from a gap so runs don't wrap
            order = [(start + k) % len(info) for k in range(len(info))]
            groups, cur = [], []
            for i in order:
                if valid[i]:
                    cur.append(i)
                elif cur:
                    groups.append(cur)
                    cur = []
            if cur:
                groups.append(cur)
            closed = False
        total = self._L[-1]
        out = []
        for g in groups:
            if len(g) < 3:
                continue
            # arc length along the run, for fading out an open run's ends
            along = [0.0]
            for a, b in zip(g, g[1:]):
                along.append(along[-1] + (info[b][0] - info[a][0]) % total)
            slices = []
            for k, i in enumerate(g):
                end = 1.0 if closed else _smoothstep(min(along[k], along[-1] - along[k]) / BLEND_END_FADE_MM)
                slices.append(self._slice(*info[i][1:], end))
            out.append((closed, slices))
        return out

    def _slice(self, x, y, nx, ny, w, reach, end_fade):
        n_up, n_lo = BLEND_SLICE_PTS
        s_in = -BLEND_INSIDE_MM
        m = max(int((reach - s_in) / 0.05), 40)
        ss = [s_in + (reach - s_in) * k / m for k in range(m + 1)]
        up, lo = [], []
        for s in ss:
            qx, qy = x + nx * s, y + ny * s
            t = self.top(qx, qy)
            f = self._fade(qx, qy, self._on_foil(qx, qy)) * end_fade
            up.append((s, t + f * self._bump(s, w) - (1.0 - f) * BLEND_SINK_MM))
            lo.append((s, max(t - BLEND_DEPTH_MM, self.low(qx, qy) + 0.3)))
        pt = lambda s, z: (x + nx * s, y + ny * s, z)
        upper = [pt(s, z) for s, z in self._even(up, n_up)]
        lower = [pt(s, z) for s, z in self._even(list(reversed(lo)), n_lo)]
        return x, y, nx, ny, upper, lower

    def rim_round_limit(self):
        """Roughly the largest round that fits on the collar's top edge. That
        edge is a wedge as sharp as the blend angle, and a round of radius r
        on it reaches r / tan(angle / 2) down each face, which has to stay
        within the collar's height."""
        return 0.8 * self.h * math.tan(math.radians(max(self.wall_deg, 1.0)) / 2.0)

    def summary(self):
        return {"blend_rise": self.h, "blend_width": self.w,
                "blend_front_width": self.w_front, "blend_back_width": self.w_back}
