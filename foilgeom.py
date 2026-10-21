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
DEFAULTS = dict(PRESETS["Bonetti main foil (E393)"])


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

    def stations(self):
        n = int(self.p["n_half"])
        eta_tip = float(self.p["tip_eta"])
        out = []
        for k in range(-n, n + 1):
            th = (math.pi / 2.0) * k / n
            out.append(math.sin(th) * eta_tip * self.b2)
        if self.bulb_on:
            yb = 1.15 * self.bulb_w
            nb = 8
            out = [y for y in out if abs(y) > yb * 1.05]
            out += [yb * k / nb for k in range(-nb, nb + 1)]
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
