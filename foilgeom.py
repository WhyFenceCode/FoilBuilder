"""foilgeom: geometry core for the FoilBuilder Fusion add-in (numpy).

Fusion's bundled Python doesn't ship numpy: FoilBuilder.py puts the add-in's
own lib/ folder (numpy installed there) on sys.path before importing this.
All lengths in mm. Point sets are numpy arrays, one point per row.
Axes: x aft (chordwise), y spanwise, z up. Origin on the planform axis at the
centreline, on the span line.

Planform: bi-elliptical. LE and TE are semi-ellipses sharing a spanwise axis
placed DELTA * c_root ahead of the root TE:
    c(y) = c_root * sqrt(1 - (2y/b)^2),    c_root = 4 A / (pi b)
Anhedral (optional, "banana" foil): span line is an elliptical arc
    z(y) = -E_V * (1 - sqrt(1 - (y/E_H)^2)),  E_H solved for the tip drop.
Sections are rotated about x to stay perpendicular to the span line.
Flap (optional): one piece behind a straight hinge line at FLAP_FRAC * c_root
ahead of the root TE, running through the bulb (so it carries the bulb's tail)
and running out each side until it tapers to nothing where the TE crosses
the hinge line.
Sections (optional): SectionedAirfoil blends several airfoils placed at
fractions of the semi-span, so the section shape can change along the span.
"""
import math

import numpy as np

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
    """Return (name, pts): pts an (n, 2) array in Selig order (TE -> upper ->
    LE -> lower -> TE), normalised to unit chord. Accepts Selig and Lednicer."""
    name = "airfoil"
    rows = []
    for i, l in enumerate(l.strip() for l in text.splitlines()):
        if not l:
            continue
        f = _floats(l)
        if f is None or len(f) < 2:
            if i == 0 or not rows:
                name = l
            continue
        rows.append(f[:2])
    if not rows:
        raise ValueError("No coordinate rows found in .dat file")
    num = np.array(rows, dtype=float)
    # Lednicer: first numeric row holds point counts (>1)
    if num[0, 0] > 1.5 and num[0, 1] > 1.5:
        n_up, n_lo = int(round(num[0, 0])), int(round(num[0, 1]))
        coords = num[1:]
        if len(coords) < n_up + n_lo:
            raise ValueError("Lednicer file: fewer points than declared")
        upper = coords[:n_up]              # LE -> TE
        lower = coords[n_up:n_up + n_lo]   # LE -> TE
        pts = np.vstack([upper[::-1], lower[1:]])
    else:
        pts = num
    # drop consecutive duplicates
    keep = np.r_[True, np.any(np.abs(np.diff(pts, axis=0)) > 1e-9, axis=1)]
    pts = pts[keep]
    xmin, xmax = pts[:, 0].min(), pts[:, 0].max()
    c = xmax - xmin
    if c <= 0:
        raise ValueError("Degenerate airfoil (zero chord)")
    pts = np.column_stack([(pts[:, 0] - xmin) / c, pts[:, 1] / c])
    # orientation check: upper surface first (Selig). If the first half is
    # below the second half, reverse.
    i_le = int(np.argmin(pts[:, 0]))
    if pts[:i_le + 1, 1].mean() < pts[i_le:, 1].mean():
        pts = pts[::-1].copy()
    # A closed TE (first point == last point) is fine: the TE opening applied in
    # Planform.section() separates the two points.
    return name, pts


def _sorted_xy(pts):
    """Rows of pts sorted by x, then y (as sorting (x, y) tuples would)."""
    return pts[np.lexsort((pts[:, 1], pts[:, 0]))]


def airfoil_stats(pts):
    """(t/c, x of max t, max camber, x of max camber) for unit-chord Selig points."""
    i_le = int(np.argmin(pts[:, 0]))
    up, lo = _sorted_xy(pts[:i_le + 1]), _sorted_xy(pts[i_le:])
    x = np.linspace(0.0, 1.0, 401)
    yu, yl = np.interp(x, up[:, 0], up[:, 1]), np.interp(x, lo[:, 0], lo[:, 1])
    out = []
    for v in (yu - yl, 0.5 * (yu + yl)):
        k = int(np.argmax(v))
        out += [float(v[k]), float(x[k])] if v[k] > 0 else [0, 0]
    return tuple(out)


# ------------------------------------------------------------------ spanwise sections

N_RESAMPLE = 50     # points per surface after resampling (cosine spaced)


def _cubic_spline(xs, ys):
    """Natural cubic spline through (xs, ys), xs strictly increasing. Returns f,
    which takes a number or an array."""
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    n = len(xs)
    if n < 3:
        return lambda x: np.interp(x, xs, ys)
    h = np.diff(xs)
    # second derivatives: tridiagonal system, natural ends (m = 0)
    A = np.zeros((n, n))
    A[0, 0] = A[-1, -1] = 1.0
    i = np.arange(1, n - 1)
    A[i, i - 1], A[i, i], A[i, i + 1] = h[:-1], 2.0 * (h[:-1] + h[1:]), h[1:]
    d = np.zeros(n)
    d[1:-1] = 6.0 * (np.diff(ys[1:]) / h[1:] - np.diff(ys[:-1]) / h[:-1])
    m = np.linalg.solve(A, d)

    def f(x):
        x = np.asarray(x, dtype=float)
        lo = np.clip(np.searchsorted(xs, x, side="right") - 1, 0, n - 2)
        x = np.clip(x, xs[lo], xs[lo + 1])
        a = (xs[lo + 1] - x) / h[lo]
        b = 1.0 - a
        return a * ys[lo] + b * ys[lo + 1] + ((a ** 3 - a) * m[lo] + (b ** 3 - b) * m[lo + 1]) * h[lo] ** 2 / 6.0
    return f


def _sqrt_spline(xs, zs):
    """Spline of z against sqrt(x - xs[0]) (smooth round a rounded leading
    edge), for x-sorted xs. Returns f(x), taking a number or an array."""
    x0 = xs[0]
    u = np.sqrt(np.maximum(xs - x0, 0.0))
    keep = np.r_[True, np.diff(u) >= 1e-9]
    f = _cubic_spline(u[keep], zs[keep])
    return lambda x: f(np.sqrt(np.maximum(np.asarray(x, dtype=float) - x0, 0.0)))


def resample(pts, n=N_RESAMPLE):
    """Resample unit-chord Selig points onto common cosine-spaced x stations, so
    that different airfoils can be blended point by point. Each surface is
    splined as y(sqrt(x)), which stays smooth round a rounded leading edge."""
    i_le = int(np.argmin(pts[:, 0]))
    xs = 0.5 * (1.0 - np.cos(np.pi * np.arange(n + 1) / n))
    up, lo = (_sorted_xy(part) for part in (pts[:i_le + 1], pts[i_le:]))
    yu = _sqrt_spline(np.maximum(up[:, 0], 0.0), up[:, 1])(xs)
    yl = _sqrt_spline(np.maximum(lo[:, 0], 0.0), lo[:, 1])(xs)
    return np.vstack([np.column_stack([xs, yu])[::-1], np.column_stack([xs, yl])[1:]])


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
        self.afs = [resample(np.asarray(s[2], dtype=float), n) for s in secs]
        self._cache = {}    # keeps blended arrays alive (Planform caches stats by id)

    @property
    def name(self):
        return " / ".join(dict.fromkeys(self.names))

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
            i = int(np.searchsorted(e, eta, side="right")) - 1
            t = (eta - e[i]) / (e[i + 1] - e[i])
            af = self.afs[i].copy()
            af[:, 1] = (1.0 - t) * self.afs[i][:, 1] + t * self.afs[i + 1][:, 1]
        self._cache[key] = af
        return af


# ------------------------------------------------------------------ parameters

TE_MIN_MM = 0.05

PRESETS = {
    "Bonetti main foil (E393)": dict(span=1150.0, area_cm2=960.0, delta=0.43, flap_on=True,
                                     flap_frac=0.30, anhedral_depth=0.0,
                                     anhedral_semiaxis=65.0, n_half=18, tip_eta=0.998,
                                     te_thick_pct=0.4, bulb_on=True, bulb_length=200.0,
                                     bulb_nose_frac=0.6, bulb_thickness=24.0, bulb_half_width=30.0),
    "Bonetti rudder foil (E393)": dict(span=800.0, area_cm2=400.0, delta=0.30, flap_on=False,
                                       flap_frac=0.30, anhedral_depth=40.0,
                                       anhedral_semiaxis=65.0, n_half=18, tip_eta=0.998,
                                       te_thick_pct=0.4, bulb_on=True, bulb_length=170.0,
                                       bulb_nose_frac=0.6, bulb_thickness=17.0, bulb_half_width=25.0),
}
STRUT_DEFAULTS = dict(strut_on=False, strut_chord=80.0, strut_le_x=25.0, strut_height=150.0,
                      strut_rim_round=0.5, strut_slot_depth=5.0, strut_flap_up_deg=15.0,
                      strut_relief_clear=1.0, strut_relief_round=5.0)
JOIN_DEFAULTS = dict(join_strut=False, join_flap=False)   # union the pieces into one body at the end
FLAP_DEFAULTS = dict(flap_angle=0.0)   # deg the flap is built deflected, + = TE down
FLAP_ANGLE_MAX = 45.0
RELIEF_GRID_MM = 0.25         # x bins of the flap's swept envelope
RELIEF_STEP_MM = 1.0          # spacing of the relief outline handed to CAD
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
BLEND_GRID_MM = 0.05          # x spacing of the foil-surface table the blend looks up
DEFAULTS = dict(PRESETS["Bonetti main foil (E393)"])
DEFAULTS.update(STRUT_DEFAULTS)
DEFAULTS.update(BLEND_DEFAULTS)
DEFAULTS.update(JOIN_DEFAULTS)
DEFAULTS.update(FLAP_DEFAULTS)


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
        if abs(float(p.get("flap_angle", 0.0))) > FLAP_ANGLE_MAX:
            raise ValueError(f"Flap angle must be within +-{FLAP_ANGLE_MAX:g} deg")
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
        """Span stations (mm), sorted array. extra_eta: semi-span fractions that
        must get a station of their own (e.g. where a user-defined section
        sits); they replace any regular station closer than 1% of the semi-span."""
        n = int(self.p["n_half"])
        eta_tip = float(self.p["tip_eta"])
        out = np.sin(0.5 * np.pi * np.arange(-n, n + 1) / n) * eta_tip * self.b2
        y_in = 0.0
        if self.bulb_on:
            yb = 1.15 * self.bulb_w
            nb = 8
            y_in = yb * 1.05
            out = np.r_[out[np.abs(out) > y_in], yb * np.arange(-nb, nb + 1) / nb]
        gap = 0.01 * self.b2
        for e in extra_eta:
            ye = abs(e) * self.b2
            if ye <= max(y_in, gap) or ye >= eta_tip * self.b2 - gap:
                continue            # already covered by centreline / bulb / tip stations
            out = np.r_[out[np.abs(np.abs(out) - ye) > gap], -ye, ye]
        return np.sort(out)

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
        """(x_hinge, y_flap_half) or None. The flap runs out along its straight
        hinge until it tapers to nothing, where the TE crosses the hinge line
        (or to the tip, for a hinge at or ahead of the ellipse axis)."""
        if not self.p.get("flap_on") or float(self.p["flap_frac"]) <= 0:
            return None
        x_te_root = self.delta * self.c_root
        x_h = x_te_root - float(self.p["flap_frac"]) * self.c_root
        s = max(x_h, 0.0) / x_te_root
        return x_h, self.b2 * math.sqrt(1.0 - s * s)

    def flap_cutter_polygon(self):
        """Planform polygon (x, y) in mm of the flap cutter, or None: everything
        aft of the hinge line, past both ends of the foil, so one flap that runs
        out as far as the hinge allows and also takes the tail of the bulb."""
        f = self.flap()
        if not f:
            return None
        x_h = f[0]
        x_far = self.delta * self.c_root + self.bulb_dL + 300.0
        y_out = self.b2 + 20.0
        return [(x_h, -y_out), (x_far, -y_out), (x_far, y_out), (x_h, y_out)]

    def flap_pivot(self, af0):
        """(x, z) of the flap's hinge axis, which runs along y: the top of the
        joint between main element and flap, on the centreline (af0: the
        centreline airfoil). The flap is rigid, so this one straight axis
        serves the whole span; away from the centreline the joint's top edge
        may sit a little lower (thinner wing, below the bulb's crown)."""
        x_h = self.flap()[0]
        (ux, uz), _ = _surfaces(self.section(af0, 0.0))
        return x_h, float(_row_spline(ux, uz)[0](x_h))

    def flap_angle(self):
        """Angle the flap is built at, deg (+ = TE down), or 0 without a flap."""
        return float(self.p.get("flap_angle", 0.0)) if self.flap_cutter_polygon() else 0.0

    def section(self, af, y):
        """(n, 3) array (mm) of the section at span station y, Selig order."""
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
        xa, za = af[:, 0], af[:, 1]
        sgn = np.sign(int(np.argmin(xa)) - np.arange(len(af)))   # +1 upper, 0 LE, -1 lower
        x_le = -(1.0 - self.delta) * c_w - d_fwd
        z0, dz = self.z_span(y)
        phi = math.atan(dz)
        zz = za * zs + sgn * xa * te / 2.0 * c
        return np.column_stack([x_le + xa * c, y - zz * math.sin(phi), z0 + zz * math.cos(phi)])

    def summary(self):
        f = self.flap()
        s = {"c_root": self.c_root, "AR": self.span ** 2 / self.area,
             "x_le_root": -(1 - self.delta) * self.c_root, "x_te_root": self.delta * self.c_root}
        if self.e_h is not None:
            y_t = float(self.p["tip_eta"]) * self.b2
            s["tip_drop"] = -self.z_span(self.b2)[0]
            s["tip_slope_deg"] = math.degrees(math.atan(abs(self.z_span(y_t)[1])))
        if f:
            s["x_hinge"] = f[0]
            s["flap_span"] = 2 * min(f[1], float(self.p["tip_eta"]) * self.b2)
            s["flap_to_tip"] = f[1] >= float(self.p["tip_eta"]) * self.b2
            s["flap_angle"] = self.flap_angle()
        if self.bulb_on:
            s["bulb_length"] = self.c_root + self.bulb_dL
            s["bulb_nose_ahead_of_le"] = self.bulb_dL * self.bulb_nose
            s["bulb_tail_x"] = self.bulb_tail_x()
            s["bulb_thickness"] = self.bulb_t
        return s


# ------------------------------------------------------------------ foil surfaces


def _surfaces(pts):
    """(n, 3) section points, Selig order -> upper and lower surface as
    x-sorted (xs, zs) arrays."""
    i_le = int(np.argmin(pts[:, 0]))
    out = []
    for part in (pts[:i_le + 1], pts[i_le:]):
        s = _sorted_xy(part[:, [0, 2]])
        out.append((s[:, 0], s[:, 1]))
    return out


def _row_spline(xs, zs):
    """Smooth z(x) along one surface of a section: (f, x_le, x_te), f taking a
    number or an array. Splined in sqrt(x - x_le), which stays smooth round a
    rounded leading edge."""
    return _sqrt_spline(xs, zs), float(xs[0]), float(xs[-1])


def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


# ------------------------------------------------------------------ vertical support


class Strut:
    """Vertical support leaving the top of the foil at the centreline.

    A straight, separate piece. Foil axes: chord along x, thickness along y,
    straight up along z to strut_height (measured from the foil chord line at
    the centreline), where it ends in a flat top face. It can sit anywhere
    along the foil, also overhanging its nose or tail.

    Its bottom is the tenon line (tenon_pts, in the x-z plane): over the foil
    it runs strut_slot_depth below the foil top, but never closer than
    STRUT_MIN_FLOOR_MM to the foil bottom (where the foil is thinner it rises
    to the surface); past the foil's nose or tail it is the height of that
    end, so the overhang is trimmed level with it instead of hanging below.
    The CAD side trims the vertical to this line and cuts it out of the foil
    as the slot, whose rim can be rounded by strut_rim_round.
    strut_slot_depth 0: no slot. The vertical sits on the foil, its bottom
    following the foil's top surface exactly (CAD trims it with the foil
    body itself); the tenon line then only trims an overhang.
    Where it reaches back over the flap, a rounded relief is cut out of its
    back so the flap can deflect up into it (see _flap_relief); the slot then
    stops ahead of the hinge, so the vertical only holds on the main element.
    eta = 0 at the highest point of the foil top under the vertical, 1 at its top."""

    def __init__(self, p, pf, foil_af_for, strut_afs, n=200):
        self.c = float(p["strut_chord"])
        self.h = float(p["strut_height"])
        self.rim = float(p.get("strut_rim_round", 0.0))
        self.slot = float(p["strut_slot_depth"])
        if self.c <= 0 or self.slot < 0 or self.rim < 0:
            raise ValueError("Vertical chord must be positive, the slot depth and rim round 0 or more")
        # strut_le_x: LE position aft of the foil nose (bulb nose with a bulb) on the centreline
        self.x_le = float(_surfaces(pf.section(foil_af_for(0.0), 0.0))[0][0][0]) + float(p["strut_le_x"])
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
        xs = np.linspace(self.x_le, self.x_te, n + 1)
        on = xs[(xs >= le0) & (xs <= te0)]
        self.z_root = float(top0(on).max())
        self.overhang = (max(0.0, le0 - self.x_le), max(0.0, self.x_te - te0))

        # tenon line, sampled ~every mm a little beyond both ends of the vertical
        m = max(int(self.c + 4.0), 8)
        xs_t = np.linspace(self.x_le - 2.0, self.x_te + 2.0, m + 1)
        self.tenon_pts = np.column_stack([xs_t, self._tenon_z(xs_t)])
        self.z_bot = float(self.tenon_pts[:, 1].min()) - 1.0       # loft bottom, below the line
        # deepest slot (on the centreline), and the least foil left under it on
        # any slice where there is a slot at all. Normally >= STRUT_MIN_FLOOR_MM;
        # less only right at a blunt nose / tail the vertical overhangs, where
        # the line has to come down to that end's height to leave no gap.
        if self.slot > 0:
            self.slot_max = float((top0(on) - self._tenon_z(on)).max())
            x, z = self.tenon_pts[:, 0], self.tenon_pts[:, 1]
            floors = [z + 0.1 - low(x) for (top, le, te), (low, _, _) in self._rows]
            mask = [(x >= le) & (x <= te) & (z < top(x) - 0.2) for (top, le, te), _ in self._rows]
            vals = np.concatenate([f[k] for f, k in zip(floors, mask)])
            self.floor = float(vals.min()) if len(vals) else 0.0
        else:
            self.slot_max = self.floor = 0.0                   # sits on the foil: no slot

        h_min = self.z_root + 5.0
        if self.h < h_min:
            raise ValueError(f"Vertical height must be at least {h_min:.0f} mm (foil top + 5 mm)")

        # flap relief: where the vertical reaches back over the flap, the part
        # of it the flap sweeps through (deflecting up) is cut away
        self.relief = None
        self.relief_x0 = None
        if pf.flap_cutter_polygon():
            self._flap_relief(p, pf.flap()[0])

    def _flap_relief(self, p, x_h):
        """Side-view (x, z) outline of the cutout in the vertical's back that
        lets the flap deflect up into it. The flap's profile under the
        vertical (highest across its thickness) is rotated about the flap
        pivot (the top of the joint on the centreline, see flap_pivot) from 0
        up to strut_flap_up_deg, or to the flap's set angle if that is further
        up; the cutout is everything under that swept
        envelope plus strut_relief_clear. Its front corner is rounded by a
        quarter circle of radius strut_relief_round that reaches forward from
        the hinge, so there is no sharp inside corner in the vertical.
        Sets self.relief ((n, 2) array, front to back; the CAD side closes it
        below self.relief_zb), self.relief_x0 (front of the cutout), or leaves
        them None when the vertical stays clear of the flap."""
        up = math.radians(float(p.get("strut_flap_up_deg", 15.0)))
        clr = float(p.get("strut_relief_clear", 1.0))
        R = float(p.get("strut_relief_round", 5.0))
        if not 0.0 <= up < math.radians(80.0) or clr < 0 or R < 0:
            raise ValueError("Flap relief: up angle 0-80 deg, clearance and round 0 or more")
        up = max(up, -math.radians(float(p.get("flap_angle", 0.0))))   # a flap set up counts too
        (top0, le0, te0), (low0, _, _) = self._rows[0]
        if not le0 < x_h < te0:
            return
        # flap silhouette under the vertical: front face, then the top surface
        xs = np.linspace(x_h, te0, max(int((te0 - x_h) / 0.2), 2) + 1)
        top = np.max([f(np.clip(xs, le, te)) for (f, le, te), _ in self._rows], axis=0)
        z_lo = float(low0(x_h))
        z_h = float(top0(x_h))                              # pivot: top of the joint, centreline
        n_f = max(int((top[0] - z_lo) / 0.2), 2)
        cx = np.r_[np.full(n_f, x_h), xs]
        cz = np.r_[np.linspace(z_lo, top[0], n_f, endpoint=False), top]
        th = np.linspace(0.0, up, max(int(math.ceil(math.degrees(up) / 0.5)), 1) + 1)[:, None]
        dx, dz = cx - x_h, cz - z_h
        X = (x_h + dx * np.cos(th) - dz * np.sin(th)).ravel()
        Z = (z_h + dx * np.sin(th) + dz * np.cos(th)).ravel()
        x_a = float(X.min()) - clr                          # front of the swept region + clearance
        if self.x_te <= x_a:
            return                                          # vertical ends ahead of the flap
        # upper envelope on an x grid (max per bin, gaps interpolated, constant
        # beyond either end), plus clearance; then made a little conservative
        # by a running max before thinning it out for CAD
        x_end = max(self.x_te, float(X.max())) + 3.0
        gx = np.arange(x_a, x_end + RELIEF_GRID_MM, RELIEF_GRID_MM)
        env = np.full(len(gx), -np.inf)
        np.maximum.at(env, np.clip(((X - x_a) / RELIEF_GRID_MM).astype(int), 0, len(gx) - 1), Z)
        ok = np.isfinite(env)
        env = np.interp(gx, gx[ok], env[ok]) + clr
        w = int(RELIEF_STEP_MM / RELIEF_GRID_MM)
        pad = np.r_[np.full(w, env[0]), env, np.full(w, env[-1])]
        env = np.max([pad[k:k + len(env)] for k in range(2 * w + 1)], axis=0)
        k = np.unique(np.r_[np.arange(0, len(gx), w), len(gx) - 1])
        line = np.column_stack([gx[k], env[k]])
        e0 = float(env[0])
        if R > 0:                                           # rounded front: vertical tangent at x_a - R
            a = np.linspace(math.pi, 0.5 * math.pi, 10)[:-1]
            arc = np.column_stack([x_a + R * np.cos(a), e0 - R + R * np.sin(a)])
            line = np.vstack([arc, line])
        self.relief_x0 = x_a - R
        self.relief_zb = min(self.z_bot, e0 - R) - 5.0
        if self.relief_x0 < self.x_le + 5.0:
            raise ValueError("The flap relief would cut away the front of the vertical: move it "
                             "forward, shorten its chord, or reduce the flap relief settings")
        if float(env.max()) > self.h - 5.0:
            raise ValueError("The flap relief reaches the vertical's top: make it taller")
        self.relief = line
        self.relief_rise = float(env[gx <= self.x_te].max()) - float(top[0])   # above the foil top at the hinge

    def _surf_at(self, x):
        """(lowest foil top, highest foil bottom) at x (array) across the
        vertical's thickness. Where a row is off the foil, its nose / tail
        height counts as both, so the line there lies level with that end."""
        tops, lows = [], []
        for (top, le, te), (low, _, _) in self._rows:
            on = (x >= le) & (x <= te)
            end = np.where(x < le, top(le), top(te))
            tops.append(np.where(on, top(x), end))
            lows.append(np.where(on, low(x), end))
        return np.min(tops, axis=0), np.max(lows, axis=0)

    def _tenon_z(self, x):
        """Height of the vertical's bottom at x: slot_depth below the foil top
        but at least STRUT_MIN_FLOOR_MM above its bottom, and never above the
        top (no gap at the vertical's sides). 0.1 mm into the foil."""
        t, b = self._surf_at(np.asarray(x, dtype=float))
        return np.minimum(t, np.maximum(t - self.slot, b + STRUT_MIN_FLOOR_MM)) - 0.1

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
        """(n, 3) array (mm) of the horizontal section at height z: Selig
        order, then a few points round a semicircular TE, as one closed loop
        (the CAD side closes it with a closed spline)."""
        xa, ya = af[:, 0], af[:, 1]
        sgn = np.sign(int(np.argmin(xa)) - np.arange(len(af)))   # +1 upper, 0 LE, -1 lower
        x = self.x_le + xa * self.c
        y = (ya + sgn * xa * self.te / 2.0) * self.c
        xm, ym = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
        r = 0.5 * math.hypot(x[0] - x[-1], y[0] - y[-1])
        t = np.radians([-60.0, -30.0, 0.0, 30.0, 60.0])          # lower TE -> aft -> upper TE
        x, y = np.r_[x, xm + r * np.cos(t)], np.r_[y, ym + r * np.sin(t)]
        return np.column_stack([x, y, np.full(len(x), float(z))])

    def summary(self):
        s = {"strut_t": 2 * self.half_t, "strut_floor": self.floor, "strut_slot_max": self.slot_max,
             "strut_z_root": self.z_root, "strut_overhang_fwd": self.overhang[0],
             "strut_overhang_aft": self.overhang[1]}
        if self.relief is not None:
            s["strut_relief_len"] = self.x_te - self.relief_x0
            s["strut_relief_rise"] = self.relief_rise
        return s


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

    N_ROWS = 25         # foil surface rows at |y| = 0 .. 24 dy (plus one each side)

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

        # the vertical's outline as a closed loop (first point repeated at the
        # end), with cumulative arc length
        P = st.section(strut_af_at(0.0), 0.0)[:, :2]
        Q = np.roll(P, -1, axis=0)
        area = np.sum(P[:, 0] * Q[:, 1] - Q[:, 0] * P[:, 1])
        self._ccw = 1.0 if area > 0 else -1.0       # (ty, -tx) * this points out
        self._loop = np.vstack([P, P[:1]])
        self._L = np.r_[0.0, np.cumsum(np.hypot(*np.diff(self._loop, axis=0).T))]

        # foil top and bottom at |y| stations round the vertical, each a smooth
        # spline in x (the blend meets the foil nearly tangent, so kinks in a
        # piecewise-linear surface would show as wobbles in its edge)
        c = np.linspace(0.0, 1.0, 51)
        side = np.max((self.w + c ** BLEND_END_POWER * (max(self.w_front, self.w_back) - self.w))
                      * np.sqrt(1.0 - c * c))
        self.y_max = st.half_t + max(side, self.w) + BLEND_SKIRT_MM + 1.0
        self._dy = (self.y_max + 1.0) / (self.N_ROWS - 1)
        self._top, self._low = [], []
        for k in range(-1, self.N_ROWS + 2):
            (ux, uz), (lx, lz) = _surfaces(pf.section(foil_af_for(k * self._dy), k * self._dy))
            self._top.append(_row_spline(ux, uz))
            self._low.append(_row_spline(lx, lz))
        self._ends = np.array([[le, te] for _, le, te in self._top])    # per row: LE x, TE x
        # every row evaluated once on a fine x grid covering anywhere a slice
        # can reach; lookups then interpolate this table in one vectorized go
        reach = max(self.w, self.w_front, self.w_back) + BLEND_SKIRT_MM + 2.0
        n = int((st.c + 2.0 * reach) / BLEND_GRID_MM) + 2
        self._gx = np.linspace(st.x_le - reach, st.x_te + reach, n)
        self._gdx = self._gx[1] - self._gx[0]
        self._T = np.array([f(self._gx) for f, _, _ in self._top])     # (rows, n)
        self._B = np.array([f(self._gx) for f, _, _ in self._low])

        # where the flap is split off, or the vertical's flap relief starts
        # (the blend fades out before it)
        self._x_stop = None
        if pf.flap_cutter_polygon():
            self._x_stop = st.relief_x0 if st.relief is not None else pf.flap()[0]

        self._stations = self._place_stations()

    # ---- foil surfaces near the centreline: rows at y = -dy .. (N_ROWS + 1) dy,
    #      each a spline in x, smooth in y by Catmull-Rom across rows
    def surface(self, x, y):
        """(top, bottom, thickness) at points (x, y) (arrays): thickness is 0
        off the foil. Linear in x on the row table, Catmull-Rom across rows."""
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        t = np.minimum(np.abs(y) / self._dy, self.N_ROWS - 1e-6) + 1.0
        k = t.astype(int)
        f = t - k
        w = np.stack([((-f + 2.0) * f - 1.0) * f, (3.0 * f - 5.0) * f * f + 2.0,
                      ((-3.0 * f + 4.0) * f + 1.0) * f, (f - 1.0) * f * f]) * 0.5   # Catmull-Rom weights
        rows = np.stack([k - 1, k, k + 1, k + 2])
        u = np.clip((x - self._gx[0]) / self._gdx, 0.0, len(self._gx) - 1.000001)
        i = u.astype(int)
        a = u - i
        pick = lambda T: np.sum(w * (T[rows, i] * (1.0 - a) + T[rows, i + 1] * a), axis=0)
        top, low = pick(self._T), pick(self._B)
        le = np.sum(w * self._ends[rows, 0], axis=0)
        te = np.sum(w * self._ends[rows, 1], axis=0)
        thick = np.where((x > le) & (x < te), top - low, 0.0)
        return top, low, thick

    def top(self, x, y):
        return self.surface(x, y)[0]

    def low(self, x, y):
        return self.surface(x, y)[1]

    # ---- the vertical's outline
    def _point(self, s):
        """Point(s) on the outline at arc length s (wraps round): x, y arrays."""
        s = np.mod(s, self._L[-1])
        return np.interp(s, self._L, self._loop[:, 0]), np.interp(s, self._L, self._loop[:, 1])

    def _normal(self, s, ds=0.1):
        """Outward unit normal(s) of the outline at arc length s (from a short
        chord, so the closely spaced points round its TE don't make it jump)."""
        (x0, y0), (x1, y1) = self._point(np.asarray(s) - ds), self._point(np.asarray(s) + ds)
        tx, ty = x1 - x0, y1 - y0
        L = np.hypot(tx, ty)
        L = np.where(L > 0, L, 1.0)
        return self._ccw * ty / L, -self._ccw * tx / L

    def _place_stations(self):
        """Arc lengths of the slices round the outline: every BLEND_SIDE_STEP_MM
        along the sides, every BLEND_FAN_STEP_DEG of turning round the ends."""
        total = self._L[-1]
        n = max(int(total / 0.05), 200)
        s = total * np.arange(n + 1) / n
        ang = np.arctan2(*self._normal(s)[::-1])
        da = np.abs((np.diff(ang) + np.pi) % (2 * np.pi) - np.pi)
        meas = np.r_[0.0, np.cumsum((total / n) / BLEND_SIDE_STEP_MM + da / math.radians(BLEND_FAN_STEP_DEG))]
        m = max(int(round(meas[-1])), 8)
        return np.interp(meas[-1] * np.arange(m) / m, meas, s)

    # ---- shape
    def _width(self, nx):
        """Blend width where the outline's outward normal has x-component nx."""
        f = np.maximum(-nx, 0.0) ** BLEND_END_POWER
        b = np.maximum(nx, 0.0) ** BLEND_END_POWER
        return self.w + f * (self.w_front - self.w) + b * (self.w_back - self.w)

    def _bump(self, d, w):
        """Height above the foil top at signed distance d (array) from the wall,
        for width w, before fading."""
        m0 = max(-w * self.cot_wall / self.h, -3.0)             # Hermite end slopes (in d / w)
        m1 = -self.tan_foot * w / self.h
        t = np.clip(d / w, 0.0, 1.0)
        collar = self.h * ((2 * t ** 3 - 3 * t ** 2 + 1) + (t ** 3 - 2 * t ** 2 + t) * m0
                           + (t ** 3 - t ** 2) * m1)
        # inside the vertical (cut away by the slot): keep rising, levelling
        # off smoothly, so there is no corner whose ripples reach the outside
        inside = self.h * (1.0 - m0 / w * (1.0 - np.exp(np.minimum(d, 0.0))))
        return np.where(d <= 0.0, inside, np.where(d >= w, -self.tan_foot * (d - w), collar))

    def _fade(self, x, y, thick):
        lo, hi = BLEND_FADE_MM
        f = _smoothstep((thick - lo) / (hi - lo))
        if self._x_stop is not None:
            f = f * _smoothstep((self._x_stop - x - 1.0) / 4.0)
        return f

    @staticmethod
    def _even(s, z, n):
        """n points evenly spaced along the polyline (s, z) (arrays): a loft
        pairs up points by their fraction along each slice, so every slice gets
        the same number, spaced the same way. Returns (s, z) arrays."""
        L = np.r_[0.0, np.cumsum(np.hypot(np.diff(s), np.diff(z)))]
        t = L[-1] * np.arange(n) / (n - 1)
        return np.interp(t, L, s), np.interp(t, L, z)

    def _reach(self, x, y, nx, ny, w):
        """How far out from the wall point (x, y) the slice goes: its width +
        skirt, or less where the foil gets too thin. None if the wall point
        itself isn't on enough foil."""
        ok = lambda s: self.surface(x + nx * s, y + ny * s)[2] >= BLEND_MIN_THICK_MM
        s_max = w + BLEND_SKIRT_MM
        s = np.r_[np.arange(0.0, s_max, 0.25), s_max]
        good = ok(s)
        if not good[0]:
            return None
        if good.all():
            return s_max
        i = int(np.argmin(good))                     # first step off enough foil
        a, b = s[i - 1], s[i]
        for _ in range(4):                           # narrow the edge down, 64 at a time
            g = np.linspace(a, b, 65)
            j = int(np.argmin(ok(g)))
            a, b = g[j - 1], g[j]
        return a

    def runs(self):
        """[(closed, [slice, ...])]: the slices to loft, in order round the
        vertical. One closed run when it goes all the way round, else open runs
        (their ends faded out). Each slice is (x, y, nx, ny, upper, lower):
        wall point, outward normal, upper points (wall outwards) and lower
        points (back to the wall), (n, 3) arrays (mm)."""
        sts = self._stations
        xs, ys = self._point(sts)
        nxs, nys = self._normal(sts)
        ws = self._width(nxs)
        reach = [self._reach(*a) for a in zip(xs, ys, nxs, nys, ws)]
        valid = np.array([r is not None and r >= 1.0 for r in reach])
        if not valid.any():
            return []
        n = len(sts)
        if valid.all():
            groups, closed = [np.arange(n)], True
        else:
            start = int(np.argmin(valid))            # walk from a gap so runs don't wrap
            order = (start + np.arange(n)) % n
            v = valid[order]
            edges = np.flatnonzero(np.diff(np.r_[0, v.astype(int), 0]))
            groups = [order[a:b] for a, b in zip(edges[::2], edges[1::2])]
            closed = False
        total = self._L[-1]
        out = []
        for g in groups:
            if len(g) < 3:
                continue
            # arc length along the run, for fading out an open run's ends
            along = np.r_[0.0, np.cumsum(np.mod(np.diff(sts[g]), total))]
            ends = np.ones(len(g)) if closed else \
                _smoothstep(np.minimum(along, along[-1] - along) / BLEND_END_FADE_MM)
            slices = [self._slice(xs[i], ys[i], nxs[i], nys[i], ws[i], reach[i], e) for i, e in zip(g, ends)]
            out.append((closed, slices))
        return out

    def _slice(self, x, y, nx, ny, w, reach, end_fade):
        n_up, n_lo = BLEND_SLICE_PTS
        s_in = -BLEND_INSIDE_MM
        m = max(int((reach - s_in) / 0.05), 40)
        s = np.linspace(s_in, reach, m + 1)
        qx, qy = x + nx * s, y + ny * s
        top, low, thick = self.surface(qx, qy)
        f = self._fade(qx, qy, thick) * end_fade
        z_up = top + f * self._bump(s, w) - (1.0 - f) * BLEND_SINK_MM
        z_lo = np.maximum(top - BLEND_DEPTH_MM, low + 0.3)
        pts = lambda s_, z_: np.column_stack([x + nx * s_, y + ny * s_, z_])
        upper = pts(*self._even(s, z_up, n_up))
        lower = pts(*self._even(s[::-1], z_lo[::-1], n_lo))
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
