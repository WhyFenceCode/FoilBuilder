"""FoilBuilder: Fusion add-in that builds a parametric T-foil (main or rudder
lifting foil) from any airfoil .dat file.

Buttons (Solid > Create panel):
  * Foil from .dat             - dialog: pick .dat, preset, parameters -> builds the foil
  * Foil from .dat by section  - same, but several .dat/.csv files placed along the
                                 semi-span; the section is blended between them
  * Rebuild Foil               - reads the foil_* user parameters and rebuilds
Both build dialogs can add a vertical support (strut): a straight separate body
that slots into a cutout in the foil/bulb.

Geometry maths lives in foilgeom.py (pure Python). This file only talks to
Fusion: construction planes, sketches, fitted splines, loft, flap split.
Fusion's internal length unit is cm; foilgeom works in mm (MM2CM below).
"""
import adsk.core
import adsk.fusion
import importlib
import json
import os
import sys
import math
import time
import traceback

_HERE = os.path.dirname(os.path.realpath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import foilgeom as G  # noqa: E402
# Fusion keeps one Python session for all add-ins, so a restarted add-in would
# otherwise get the foilgeom it imported before any edits. Always load it fresh.
G = importlib.reload(G)

APP_NAME = "FoilBuilder"
CMD_BUILD = "aquafoil_FoilBuilder_build"
CMD_BUILD_SEC = "aquafoil_FoilBuilder_build_sections"
CMD_REBUILD = "aquafoil_FoilBuilder_rebuild"
PANEL_ID = "SolidCreatePanel"
WORKSPACE_ID = "FusionSolidEnvironment"
ATTR_GROUP = "aQuaFoil_FoilBuilder"
MM2CM = 0.1
DEFAULT_UP_AXIS = "Y"      # used when a design has no stored choice

# user parameter name, key in foilgeom params, unit ('' = unitless), comment
PARAMS = [
    ("foil_span", "span", "mm", "Projected span"),
    ("foil_area_cm2", "area_cm2", "", "Planform area [cm^2]"),
    ("foil_delta", "delta", "", "Ellipse axis position from TE, fraction of root chord"),
    ("foil_flap_on", "flap_on", "", "1 = split a flap, 0 = no flap"),
    ("foil_flap_frac", "flap_frac", "", "Flap chord / root chord (straight hinge)"),
    ("foil_flap_min_chord", "flap_min_chord", "mm", "Flap ends where its chord drops below this"),
    ("foil_anhedral_depth", "anhedral_depth", "mm", "Tip drop of curved anhedral (0 = flat)"),
    ("foil_anhedral_semiaxis", "anhedral_semiaxis", "mm", "Vertical semi-axis of anhedral ellipse"),
    ("foil_n_half", "n_half", "", "Sections per half span"),
    ("foil_tip_eta", "tip_eta", "", "Last section position, fraction of semi-span"),
    ("foil_te_thick_pct", "te_thick_pct", "", "Trailing-edge thickness, % of local chord"),
    ("foil_bulb_on", "bulb_on", "", "1 = blended bulb at the centreline, 0 = none"),
    ("foil_bulb_length", "bulb_length", "mm", "Bulb chord at the centreline"),
    ("foil_bulb_nose_frac", "bulb_nose_frac", "", "Share of extra bulb length ahead of the wing LE"),
    ("foil_bulb_thickness", "bulb_thickness", "mm", "Bulb max thickness at the centreline"),
    ("foil_bulb_half_width", "bulb_half_width", "mm", "Spanwise half-width over which the bulb blends out"),
    ("foil_strut_on", "strut_on", "", "1 = vertical support on top of the foil, 0 = none"),
    ("foil_strut_chord", "strut_chord", "mm", "Vertical support chord"),
    ("foil_strut_le_x", "strut_le_x", "mm", "Vertical LE aft of the foil nose on the centreline"),
    ("foil_strut_height", "strut_height", "mm", "Height of the vertical's flat top above the foil chord line"),
    ("foil_strut_rim_round", "strut_rim_round", "mm", "Round on the slot rim where the vertical enters the foil (0 = sharp)"),
    ("foil_strut_slot_depth", "strut_slot_depth", "mm", "Slot depth below the lowest foil-top point under the vertical"),
]
BOOL_KEYS = ("flap_on", "bulb_on", "strut_on")

_app = None
_ui = None
_handlers = []


# ============================================================ helpers
def _design():
    d = adsk.fusion.Design.cast(_app.activeProduct)
    if not d:
        raise RuntimeError("Open a Fusion design (Design workspace) first.")
    return d


def _write_params(design, p, foil_src, up_axis, strut_src=None):
    """foil_src / strut_src: a .dat path, or [(eta, path)] for a by-section
    build. A by-section build leaves the stored single .dat paths as they were;
    a single-.dat build clears the stored sections."""
    ups = design.userParameters
    for name, key, unit, comment in PARAMS:
        val = p.get(key, G.DEFAULTS.get(key, 0))
        if key in BOOL_KEYS:
            val = 1 if val else 0
        expr = f"{val} {unit}".strip()
        existing = ups.itemByName(name)
        if existing:
            existing.expression = expr
        else:
            ups.add(name, adsk.core.ValueInput.createByString(expr), unit, comment)
    for src, single, multi in ((foil_src, "dat_path", "sections"),
                               (strut_src, "strut_dat_path", "strut_sections")):
        if isinstance(src, (list, tuple)):
            design.attributes.add(ATTR_GROUP, multi, json.dumps([[e, s] for e, s in src]))
        elif src:
            design.attributes.add(ATTR_GROUP, single, src)
    if not isinstance(foil_src, (list, tuple)):
        old = design.attributes.itemByName(ATTR_GROUP, "sections")
        if old:
            old.deleteMe()
    design.attributes.add(ATTR_GROUP, "up_axis", up_axis)


def _read_attr(design, key):
    a = design.attributes.itemByName(ATTR_GROUP, key)
    return a.value if a else ""


def _read_sections(design, key="sections"):
    """[(eta, path)] stored by the last by-section build, or [] (single .dat)."""
    v = _read_attr(design, key)
    if not v:
        return []
    try:
        return [(float(e), str(s)) for e, s in json.loads(v)]
    except Exception:
        return []


def _read_params(design):
    ups = design.userParameters
    p = dict(G.DEFAULTS)
    for name, key, unit, _ in PARAMS:
        prm = ups.itemByName(name)
        if not prm:
            continue
        v = prm.value
        if unit == "mm":
            v = v / MM2CM            # internal cm -> mm
        if key in BOOL_KEYS:
            v = v > 0.5
        elif key == "n_half":
            v = int(round(v))
        p[key] = v
    if not ups.itemByName("foil_bulb_on"):
        p["bulb_on"] = False          # designs saved before the bulb option existed
    a = design.attributes.itemByName(ATTR_GROUP, "dat_path")
    u = design.attributes.itemByName(ATTR_GROUP, "up_axis")
    return p, (a.value if a else ""), (u.value if u else DEFAULT_UP_AXIS)


def _mapper(up_axis):
    """Foil axes (x aft, y span, z up) in mm -> Fusion model Point3D in cm."""
    if up_axis == "Y":
        return lambda x, y, z: adsk.core.Point3D.create(x * MM2CM, z * MM2CM, -y * MM2CM)
    return lambda x, y, z: adsk.core.Point3D.create(x * MM2CM, y * MM2CM, z * MM2CM)


def _plane_through(comp, helper, p0, p1, p2):
    """Construction plane through three model-space points (via helper sketch points)."""
    pts = [helper.sketchPoints.add(helper.modelToSketchSpace(p)) for p in (p0, p1, p2)]
    pin = comp.constructionPlanes.createInput()
    if not pin.setByThreePoints(pts[0], pts[1], pts[2]):
        raise RuntimeError("Could not define a construction plane")
    pl = _tag(comp.constructionPlanes.add(pin))
    pl.isLightBulbOn = False
    return pl


# ---- tracking of everything a build creates (so Rebuild can remove it again,
#      also in single-component Part Design documents)
_build = {"id": "0", "n": 0}


def _tag(obj):
    """Mark an entity as created by FoilBuilder, with its creation order."""
    try:
        obj.attributes.add(ATTR_GROUP, "seq", f"{_build['id']}:{_build['n']:06d}")
        _build["n"] += 1
    except Exception:
        pass
    return obj


def _delete_tagged(design, only_build=None):
    items = []
    try:
        attrs = design.findAttributes(ATTR_GROUP, "seq")
    except Exception:
        attrs = []
    for a in attrs:
        try:
            bid, n = a.value.split(":")
            if only_build is None or bid == only_build:
                items.append((bid, n, a.parent))
        except Exception:
            pass
    for _, _, ent in sorted(items, key=lambda t: (t[0], t[1]), reverse=True):
        try:
            if ent is not None and getattr(ent, "isValid", True):
                ent.deleteMe()
        except Exception:
            pass


def _delete_previous(design):
    root = design.rootComponent
    # assembly-style builds: the foil lives in its own component
    for i in range(root.occurrences.count - 1, -1, -1):
        occ = root.occurrences.item(i)
        if occ.component.attributes.itemByName(ATTR_GROUP, "generated"):
            occ.deleteMe()
    # part-style builds: features live in the root component
    try:
        groups = design.timeline.timelineGroups
        for i in range(groups.count - 1, -1, -1):
            g = groups.item(i)
            if g.name.startswith("FoilBuilder"):
                g.deleteMe(False)          # ungroup, contents removed below
    except Exception:
        pass
    _delete_tagged(design)


# ============================================================ the builder
def _load_dat(path):
    with open(path, "r", errors="ignore") as f:
        return G.parse_dat(f.read())


def _load_sections(sections):
    """[(eta, path)] -> G.SectionedAirfoil."""
    secs = []
    for eta, path in sections:
        if not path or not os.path.isfile(path):
            raise ValueError(f"Section file not found: {path}")
        name, pts = _load_dat(path)
        secs.append((float(eta), name, pts))
    return G.SectionedAirfoil(secs)


def _airfoil_source(src):
    """src: a .dat/.csv path, or [(eta, path)] for by-section ->
    (name, af_at(eta), etas, afs). eta is the caller's span fraction."""
    if isinstance(src, (list, tuple)):
        saf = _load_sections(src)
        return saf.name, saf.at, saf.etas, saf.afs
    if not src or not os.path.isfile(src):
        raise ValueError(f"Airfoil file not found: {src}")
    name, af = _load_dat(src)
    return name, (lambda eta: af), (), [af]


def _prepare(p, foil_src, strut_src=None):
    """Load the airfoils and run all geometry checks.
    Returns (name, pf, af_for(y), extra_eta, strut) with strut = None or
    (G.Strut, af_at(eta), etas)."""
    pf = G.Planform(p)
    name, af_at, extra_eta, _ = _airfoil_source(foil_src)
    af_for = lambda y: af_at(abs(y) / pf.b2)
    strut = None
    if p.get("strut_on"):
        if not strut_src:
            raise ValueError("Pick an airfoil file for the vertical")
        _, s_at, s_etas, s_afs = _airfoil_source(strut_src)
        strut = (G.Strut(p, pf, af_for, s_afs), s_at, s_etas)
    return name, pf, af_for, extra_eta, strut


def build_foil(p, foil_src, up_axis=DEFAULT_UP_AXIS, replace=True, strut_src=None):
    """foil_src / strut_src: a .dat/.csv path, or [(eta, path)] to build by section."""
    design = _design()
    af_name, pf, af_for, extra_eta, strut = _prepare(p, foil_src, strut_src)
    P = _mapper(up_axis)
    root = design.rootComponent
    if replace:
        _delete_previous(design)

    label = f"Foil {af_name} b{p['span']:.0f} A{p['area_cm2']:.0f}"
    _build["id"] = f"{int(time.time() * 1000):015d}"
    _build["n"] = 0
    try:
        tl_start = design.timeline.markerPosition
    except Exception:
        tl_start = None

    # Assembly / classic designs: put the foil in its own component.
    # Part Design documents allow only one component: build in the root.
    try:
        occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
        comp = occ.component
        comp.name = label
        comp.attributes.add(ATTR_GROUP, "generated", "1")
    except Exception:
        occ, comp = None, root

    try:
        notes = _build_geometry(design, comp, af_for, pf, P, extra_eta, strut)
    except Exception:
        if occ is not None:
            try:
                occ.deleteMe()
            except Exception:
                pass
        else:
            _delete_tagged(design, only_build=_build["id"])   # clean up partial build
        raise

    # tidy timeline
    try:
        tl = design.timeline
        if tl_start is not None and tl.markerPosition - 1 > tl_start:
            grp = tl.timelineGroups.add(tl_start, tl.markerPosition - 1)
            grp.name = "FoilBuilder: " + label
    except Exception:
        pass

    _write_params(design, p, foil_src, up_axis, strut_src)
    s = pf.summary()
    if strut:
        s.update(strut[0].summary())
    s.update(notes)
    return af_name, s


def _section_profile(comp, plane, pts, P, name, smooth_te=False):
    """Sketch on plane: fitted spline through pts (foil axes, mm) closed by a
    TE line, or with smooth_te by closing the spline itself (one smooth loop
    with a tiny rounded TE). Returns its single profile."""
    sk = _tag(comp.sketches.add(plane))
    sk.name = name
    sk.isComputeDeferred = True
    coll = adsk.core.ObjectCollection.create()
    for (x, y, z) in pts:
        coll.add(sk.modelToSketchSpace(P(x, y, z)))
    spl = sk.sketchCurves.sketchFittedSplines.add(coll)
    if smooth_te:
        spl.isClosed = True
    else:
        sk.sketchCurves.sketchLines.addByTwoPoints(spl.endSketchPoint, spl.startSketchPoint)
    sk.isComputeDeferred = False
    if sk.profiles.count != 1:
        raise RuntimeError(f"{name} gave {sk.profiles.count} profiles")
    sk.isVisible = False
    adsk.doEvents()
    return sk.profiles.item(0)


def _loft(comp, profiles):
    lofts = comp.features.loftFeatures
    li = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    for pr in profiles:
        li.loftSections.add(pr)
    li.isSolid = True
    return _tag(lofts.add(li)).bodies.item(0)


def _combine(comp, target, tool, op, keep_tool, feature=False):
    """Combine tool into target; returns the resulting target body (or the feature)."""
    tools = adsk.core.ObjectCollection.create()
    tools.add(tool)
    ci = comp.features.combineFeatures.createInput(target, tools)
    ci.operation = op
    ci.isKeepToolBodies = keep_tool
    feat = _tag(comp.features.combineFeatures.add(ci))
    return feat if feature else feat.bodies.item(0)


def _build_strut(comp, helper, foil, strut, P):
    """Vertical support: a straight loft from the slot bottom to its flat top,
    kept as its own body. The same shape is cut out of the foil as the slot,
    and the slot rim on the foil is optionally rounded.
    Returns (slotted foil, vertical, notes)."""
    st, af_at, etas = strut
    # One smooth closed loop per section (tiny rounded TE): a flat TE strip
    # would leave two sharp corners in the slot rim, which can't be rounded.
    profiles = []
    for z in st.stations(etas):
        plane = _plane_through(comp, helper, P(0.0, 0.0, z), P(10.0, 0.0, z), P(0.0, 10.0, z))
        profiles.append(_section_profile(comp, plane, st.section(af_at(st.eta(z)), z), P,
                                         f"vertical z={z:.1f}", smooth_te=True))
    vert = _loft(comp, profiles)
    vert.name = "vertical"
    foil = _combine(comp, foil, vert, adsk.fusion.FeatureOperations.CutFeatureOperation, True)
    foil.name = "foil"

    notes = {}
    if st.rim > 0:
        try:
            _round_slot_rim(comp, foil, st, P)
        except Exception as e:
            notes["rim_round_failed"] = f"{e}".splitlines()[0]
    return foil, vert, notes


def _round_slot_rim(comp, foil, st, P):
    """Round the slot rim: edges between a slot wall (a face inside the
    vertical's footprint) and the foil's outer surface."""
    m = 1.0                                   # mm slack (face bounding boxes can be loose)
    corners = [P(x, y, z) for x in (st.x_le - m, st.x_te + m) for y in (-st.half_t - m, st.half_t + m)
               for z in (st.z_bot - m, st.h + m)]
    lo = [min(getattr(c, a) for c in corners) for a in "xyz"]
    hi = [max(getattr(c, a) for c in corners) for a in "xyz"]

    def inside(face):
        bb = face.boundingBox
        return all(lo[i] <= getattr(bb.minPoint, a) and getattr(bb.maxPoint, a) <= hi[i]
                   for i, a in enumerate("xyz"))

    edges = adsk.core.ObjectCollection.create()
    for i in range(foil.edges.count):
        e = foil.edges.item(i)
        if e.faces.count == 2 and inside(e.faces.item(0)) != inside(e.faces.item(1)):
            edges.add(e)
    if edges.count == 0:
        raise RuntimeError("slot rim edges not found")
    fillets = comp.features.filletFeatures
    fi = fillets.createInput()
    radius = adsk.core.ValueInput.createByReal(st.rim * MM2CM)
    try:
        fi.edgeSetInputs.addConstantRadiusEdgeSet(edges, radius, True)
    except AttributeError:                       # older API
        fi.addConstantRadiusEdgeSet(edges, radius, True)
    _tag(fillets.add(fi))


def _build_geometry(design, comp, af_for, pf, P, extra_eta=(), strut=None):
    """af_for(y) -> unit-chord Selig points of the airfoil at span station y.
    strut: None or (G.Strut, af_at(eta), etas) for the vertical support."""
    helper = _tag(comp.sketches.add(comp.xYConstructionPlane))
    helper.name = "FoilBuilder helper points"

    profiles = []
    for y in pf.stations(extra_eta):
        z0, dz = pf.z_span(y)
        phi = math.atan(dz)
        ny, nz = -math.sin(phi), math.cos(phi)
        plane = _plane_through(comp, helper, P(0.0, y, z0), P(10.0, y, z0),
                               P(0.0, y + 10.0 * ny, z0 + 10.0 * nz))
        profiles.append(_section_profile(comp, plane, pf.section(af_for(y), y), P,
                                         f"section y={y:.1f}"))
    foil = _loft(comp, profiles)
    foil.name = "foil"

    # ---- vertical support (before the flap split, which only touches the foil)
    notes = {}
    if strut:
        foil, _, notes = _build_strut(comp, helper, foil, strut, P)
    helper.isVisible = False

    # ---- flap split (straight hinge line)
    poly = pf.flap_cutter_polygon()
    if poly:
        x_h = pf.flap()[0]
        pl = _plane_through(comp, helper, P(0, 0, 0), P(10, 0, 0), P(0, 10, 0))
        sk = _tag(comp.sketches.add(pl))
        sk.name = "flap cutter"
        c = [sk.modelToSketchSpace(P(x, y, 0)) for (x, y) in poly]
        lines = sk.sketchCurves.sketchLines
        n = len(c)
        for i in range(n):
            lines.addByTwoPoints(c[i], c[(i + 1) % n])
        ext = comp.features.extrudeFeatures
        ei = ext.createInput(sk.profiles.item(0), adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        ei.setSymmetricExtent(adsk.core.ValueInput.createByReal(20.0), True)   # 200 mm total
        cutter = _tag(ext.add(ei)).bodies.item(0)
        sk.isVisible = False
        si = comp.features.splitBodyFeatures.createInput(foil, cutter, False)
        split = _tag(comp.features.splitBodyFeatures.add(si))
        _tag(comp.features.removeFeatures.add(cutter))
        # name the pieces: flaps start at the hinge line; with a bulb there are
        # two flap halves, named by which side of the centreline they are on
        x_h_model = x_h * MM2CM
        o, e = P(0, 0, 0), P(0, 10, 0)
        d = (e.x - o.x, e.y - o.y, e.z - o.z)                 # foil +y direction
        for i in range(split.bodies.count):
            b = split.bodies.item(i)
            bb = b.boundingBox
            if bb.minPoint.x < x_h_model - 1e-3:
                b.name = "main element"
            elif pf.bulb_on:
                cy = sum(0.5 * (a + c_) * k for a, c_, k in zip(
                    (bb.minPoint.x, bb.minPoint.y, bb.minPoint.z),
                    (bb.maxPoint.x, bb.maxPoint.y, bb.maxPoint.z), d))
                b.name = "flap starboard" if cy > 0 else "flap port"
            else:
                b.name = "flap"
    return notes


# ============================================================ dialog
AIRFOIL_FILTER = "Airfoil (*.dat;*.csv;*.txt);;All files (*.*)"

# vertical support spinners: key, label, min, max, step
STRUT_SPINNERS = (
    ("strut_chord", "Chord (mm)", 10, 400, 1),
    ("strut_le_x", "LE aft of foil nose (mm)", -100, 400, 1),
    ("strut_height", "Top height above chord line (mm)", 10, 2000, 5),
    ("strut_rim_round", "Slot rim round (mm, 0 = sharp)", 0, 10, 0.25),
    ("strut_slot_depth", "Slot depth into foil (mm)", 0.5, 50, 0.5),
)
# inputs shown only while the vertical support is ticked
STRUT_IDS = ("strut_help", "strut_table", "strut_add", "strut_remove", "strut_browse",
             "strut_dat_path") + tuple(s[0] for s in STRUT_SPINNERS)
# section table buttons -> table
ADD_BUTTONS = {"sec_add": "sec_table", "strut_add": "strut_table"}
REMOVE_BUTTONS = {"sec_remove": "sec_table", "strut_remove": "strut_table"}


def _af_line(name, af):
    t, xt, cm, xc = G.airfoil_stats(af)
    return (f"<b>{name}</b>: {len(af)} pts, t/c {t*100:.1f}% @ {xt*100:.0f}%, "
            f"camber {cm*100:.1f}% @ {xc*100:.0f}%")


def _src_lines(src):
    """Info lines for a .dat path or [(eta, path)]."""
    if isinstance(src, (list, tuple)):
        return [f"η {eta:.3f} " + _af_line(*_load_dat(path))
                for eta, path in sorted(src) if path and os.path.isfile(path)]
    return [_af_line(*_load_dat(src))] if src and os.path.isfile(src) else []


def _sources_from_inputs(inputs, by_section):
    """(foil_src, strut_src) as build_foil takes them."""
    if by_section:
        return _sections_from_inputs(inputs, "sec_table"), _sections_from_inputs(inputs, "strut_table")
    return inputs.itemById("dat_path").text, inputs.itemById("strut_dat_path").text


def _check(inputs, by_section):
    """Run every check a build would. Returns (info lines, error message or None)."""
    lines = []
    try:
        p = _params_from_inputs(inputs)
        foil_src, strut_src = _sources_from_inputs(inputs, by_section)
        lines += _src_lines(foil_src)
        s = G.Planform(p).summary()
        txt = f"Root chord {s['c_root']:.1f} mm, AR {s['AR']:.1f}"
        if "x_hinge" in s:
            txt += f", flap span {s['flap_span']:.0f} mm"
        if "tip_drop" in s:
            txt += f", tip slope {s['tip_slope_deg']:.1f} deg"
        if "bulb_length" in s:
            txt += (f"<br>Bulb {s['bulb_length']:.0f} mm long, nose {s['bulb_nose_ahead_of_le']:.0f} mm "
                    f"ahead of wing LE, {s['bulb_thickness']:.0f} mm thick")
        lines.append(txt)
        if p["strut_on"]:
            lines += ["Vertical " + l for l in _src_lines(strut_src)]
        if not foil_src:
            raise ValueError("Add at least one section" if by_section else "Pick an airfoil file")
        strut = _prepare(p, foil_src, strut_src)[4]
        if strut:
            ss = strut[0].summary()
            lines.append(f"Vertical {ss['strut_t']:.1f} mm thick, slot bottom at {ss['strut_z_bot']:.1f} mm, "
                         f"{ss['strut_floor']:.1f} mm of foil left under the slot")
        return lines, None
    except Exception as e:
        return lines, str(e)


def _info_text(inputs, by_section=False):
    lines, err = _check(inputs, by_section)
    if err:
        lines.append(f"<i>{err}</i>")
    return "<br>".join(lines)


def _params_from_inputs(inputs):
    g = lambda i: inputs.itemById(i)
    p = dict(
        span=g("span").value, area_cm2=g("area_cm2").value, delta=g("delta").value,
        flap_on=g("flap_on").value, flap_frac=g("flap_frac").value,
        flap_min_chord=g("flap_min_chord").value, anhedral_depth=g("anhedral_depth").value,
        anhedral_semiaxis=g("anhedral_semiaxis").value, n_half=g("n_half").value,
        tip_eta=g("tip_eta").value, te_thick_pct=g("te_thick_pct").value,
        bulb_on=g("bulb_on").value, bulb_length=g("bulb_length").value,
        bulb_nose_frac=g("bulb_nose_frac").value, bulb_thickness=g("bulb_thickness").value,
        bulb_half_width=g("bulb_half_width").value, strut_on=g("strut_on").value)
    for key, *_ in STRUT_SPINNERS:
        p[key] = g(key).value
    return p


def _set_inputs(inputs, p):
    for key in ("span", "area_cm2", "delta", "flap_frac", "flap_min_chord", "anhedral_depth",
                "anhedral_semiaxis", "tip_eta", "te_thick_pct", "bulb_length", "bulb_nose_frac",
                "bulb_thickness", "bulb_half_width"):
        inputs.itemById(key).value = float(p.get(key, G.DEFAULTS[key]))
    inputs.itemById("n_half").value = int(p["n_half"])
    inputs.itemById("flap_on").value = bool(p["flap_on"])
    inputs.itemById("bulb_on").value = bool(p.get("bulb_on", False))


def _show_strut(inputs):
    on = inputs.itemById("strut_on").value
    for cid in STRUT_IDS:
        it = inputs.itemById(cid)
        if it:
            it.isVisible = on


# ---- section tables (by-section dialog): one row per airfoil, [eta | file].
#      Full paths are kept here, keyed by the row's file text box id.
_sec_paths = {}
_sec_seq = [0]


def _button(inputs, cid, name, text):
    btn = inputs.addBoolValueInput(cid, name, False, "", False)
    try:
        btn.text = text
    except Exception:
        pass
    return btn


def _add_table(inputs, prefix, table_id, label, help_text):
    inputs.addTextBoxCommandInput(f"{prefix}_help", "", help_text, 2, True)
    tbl = inputs.addTableCommandInput(table_id, label, 2, "1:3")
    tbl.minimumVisibleRows = 3
    tbl.maximumVisibleRows = 8
    _button(inputs, f"{prefix}_add", "Add airfoils", "Add...")
    _button(inputs, f"{prefix}_remove", "Remove row", "Remove selected (or last)")


def _add_section_row(inputs, eta, path, table_id="sec_table"):
    tbl = inputs.itemById(table_id)
    ci = tbl.commandInputs
    k = _sec_seq[0]
    _sec_seq[0] += 1
    sp = ci.addFloatSpinnerCommandInput(f"sec_eta_{k}", "η", "", 0, 1, 0.05, eta)
    tx = ci.addTextBoxCommandInput(f"sec_file_{k}", "File", os.path.basename(path), 1, True)
    tx.tooltip = path
    _sec_paths[tx.id] = path
    row = tbl.rowCount
    tbl.addCommandInput(sp, row, 0)
    tbl.addCommandInput(tx, row, 1)


def _sections_from_inputs(inputs, table_id="sec_table"):
    tbl = inputs.itemById(table_id)
    out = []
    for r in range(tbl.rowCount):
        sp, tx = tbl.getInputAtPosition(r, 0), tbl.getInputAtPosition(r, 1)
        if sp and tx:
            out.append((sp.value, _sec_paths.get(tx.id, "")))
    return out


class BuildCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def __init__(self, by_section=False):
        super().__init__()
        self.by_section = by_section

    def notify(self, args):
        try:
            cmd = args.command
            inputs = cmd.commandInputs
            try:
                design = _design()
                p, dat, up = _read_params(design)
                secs = _read_sections(design)
                strut_dat = _read_attr(design, "strut_dat_path")
                strut_secs = _read_sections(design, "strut_sections")
            except Exception:
                p, dat, up, secs = dict(G.DEFAULTS), "", DEFAULT_UP_AXIS, []
                strut_dat, strut_secs = "", []

            if self.by_section:
                _sec_paths.clear()
                _add_table(inputs, "sec", "sec_table", "Sections",
                           "One airfoil (.dat or .csv) per row. η = position along the semi-span "
                           "(0 = centreline, 1 = tip), mirrored to both sides. The section is "
                           "blended linearly between rows.")
                for eta, path in (secs or ([(0.0, dat)] if dat else [])):
                    _add_section_row(inputs, eta, path)
            else:
                _button(inputs, "browse", "Airfoil .dat", "Browse...")
                inputs.addTextBoxCommandInput("dat_path", "File", dat, 2, True)
            dd = inputs.addDropDownCommandInput("preset", "Preset",
                                                adsk.core.DropDownStyles.TextListDropDownStyle)
            dd.listItems.add("(keep current values)", True, "")
            for k in G.PRESETS:
                dd.listItems.add(k, False, "")

            g = inputs.addGroupCommandInput("g_plan", "Planform").children
            g.addFloatSpinnerCommandInput("span", "Span (mm)", "", 100, 3000, 10, p["span"])
            g.addFloatSpinnerCommandInput("area_cm2", "Area (cm²)", "", 50, 3000, 10, p["area_cm2"])
            g.addFloatSpinnerCommandInput("delta", "Δ axis from TE (0-1)", "", 0.05, 0.95, 0.01, p["delta"])

            g = inputs.addGroupCommandInput("g_flap", "Flap").children
            g.addBoolValueInput("flap_on", "Split flap", True, "", bool(p["flap_on"]))
            g.addFloatSpinnerCommandInput("flap_frac", "Flap / root chord", "", 0.05, 0.6, 0.01, p["flap_frac"])
            g.addFloatSpinnerCommandInput("flap_min_chord", "Flap min chord (mm)", "", 1, 50, 1, p["flap_min_chord"])

            g = inputs.addGroupCommandInput("g_anh", "Curved anhedral").children
            g.addFloatSpinnerCommandInput("anhedral_depth", "Tip drop (mm, 0 = flat)", "", 0, 300, 1, p["anhedral_depth"])
            g.addFloatSpinnerCommandInput("anhedral_semiaxis", "Ellipse semi-axis (mm)", "", 1, 1000, 1, p["anhedral_semiaxis"])

            g = inputs.addGroupCommandInput("g_bulb", "Bulb (blended)").children
            g.addBoolValueInput("bulb_on", "Blended bulb", True, "", bool(p.get("bulb_on", False)))
            g.addFloatSpinnerCommandInput("bulb_length", "Length at centreline (mm)", "", 10, 600, 5,
                                          p.get("bulb_length", G.DEFAULTS["bulb_length"]))
            g.addFloatSpinnerCommandInput("bulb_nose_frac", "Share of extra length ahead of LE", "", 0, 1, 0.05,
                                          p.get("bulb_nose_frac", G.DEFAULTS["bulb_nose_frac"]))
            g.addFloatSpinnerCommandInput("bulb_thickness", "Max thickness (mm)", "", 2, 100, 0.5,
                                          p.get("bulb_thickness", G.DEFAULTS["bulb_thickness"]))
            g.addFloatSpinnerCommandInput("bulb_half_width", "Blend half-width (mm)", "", 2, 200, 1,
                                          p.get("bulb_half_width", G.DEFAULTS["bulb_half_width"]))

            g = inputs.addGroupCommandInput("g_strut", "Vertical support").children
            g.addBoolValueInput("strut_on", "Vertical support", True, "", bool(p.get("strut_on", False)))
            if self.by_section:
                _add_table(g, "strut", "strut_table", "Vertical sections",
                           "One airfoil per row. η = height from the foil top (0) to the "
                           "vertical's flat top (1).")
                for eta, path in (strut_secs or ([(0.0, strut_dat)] if strut_dat else [])):
                    _add_section_row(inputs, eta, path, "strut_table")
            else:
                _button(g, "strut_browse", "Vertical .dat", "Browse...")
                g.addTextBoxCommandInput("strut_dat_path", "File", strut_dat, 2, True)
            for key, label, lo, hi, step in STRUT_SPINNERS:
                g.addFloatSpinnerCommandInput(key, label, "", lo, hi, step,
                                              float(p.get(key, G.DEFAULTS[key])))
            _show_strut(inputs)

            g = inputs.addGroupCommandInput("g_adv", "Advanced").children
            g.addIntegerSpinnerCommandInput("n_half", "Sections per half", 4, 60, 1, int(p["n_half"]))
            g.addFloatSpinnerCommandInput("tip_eta", "Tip cut (fraction of semi-span)", "", 0.9, 0.9999, 0.001, p["tip_eta"])
            g.addFloatSpinnerCommandInput("te_thick_pct", "TE thickness (% chord)", "", 0.05, 3, 0.05, p["te_thick_pct"])
            ud = g.addDropDownCommandInput("up_axis", "Vertical axis", adsk.core.DropDownStyles.TextListDropDownStyle)
            ud.listItems.add("Y", up != "Z", "")
            ud.listItems.add("Z", up == "Z", "")
            g.addBoolValueInput("replace", "Replace previous foil", True, "", True)

            inputs.addTextBoxCommandInput("info", "", _info_text(inputs, self.by_section), 4, True)

            h1 = BuildExecuteHandler(self.by_section); cmd.execute.add(h1); _handlers.append(h1)
            h2 = BuildChangedHandler(self.by_section); cmd.inputChanged.add(h2); _handlers.append(h2)
            h3 = BuildValidateHandler(self.by_section); cmd.validateInputs.add(h3); _handlers.append(h3)
        except Exception:
            _ui.messageBox("FoilBuilder dialog failed:\n" + traceback.format_exc())


class BuildChangedHandler(adsk.core.InputChangedEventHandler):
    def __init__(self, by_section=False):
        super().__init__()
        self.by_section = by_section

    def notify(self, args):
        try:
            inputs = args.firingEvent.sender.commandInputs   # top level (args.inputs is per group)
            changed = args.input
            if changed.id in ("browse", "strut_browse"):
                dlg = _ui.createFileDialog()
                dlg.title = "Select airfoil coordinate file"
                dlg.filter = AIRFOIL_FILTER
                if dlg.showOpen() == adsk.core.DialogResults.DialogOK:
                    target = "dat_path" if changed.id == "browse" else "strut_dat_path"
                    inputs.itemById(target).text = dlg.filename
            elif changed.id in ADD_BUTTONS:
                table_id = ADD_BUTTONS[changed.id]
                dlg = _ui.createFileDialog()
                dlg.title = "Select airfoil coordinate file(s)"
                dlg.filter = AIRFOIL_FILTER
                dlg.isMultiSelectEnabled = True
                if dlg.showOpen() == adsk.core.DialogResults.DialogOK:
                    files = list(dlg.filenames)
                    n0 = inputs.itemById(table_id).rowCount
                    for i, path in enumerate(files):
                        if n0 == 0:      # fresh table: spread 0 -> 1
                            eta = i / (len(files) - 1) if len(files) > 1 else 0.0
                        else:
                            eta = 1.0
                        _add_section_row(inputs, eta, path, table_id)
            elif changed.id in REMOVE_BUTTONS:
                tbl = inputs.itemById(REMOVE_BUTTONS[changed.id])
                r = tbl.selectedRow
                if r < 0 or r >= tbl.rowCount:
                    r = tbl.rowCount - 1
                if r >= 0:
                    _sec_paths.pop(tbl.getInputAtPosition(r, 1).id, None)
                    tbl.deleteRow(r)
            elif changed.id == "strut_on":
                _show_strut(inputs)
            elif changed.id == "preset":
                name = changed.selectedItem.name if changed.selectedItem else ""
                if name in G.PRESETS:
                    _set_inputs(inputs, G.PRESETS[name])
            inputs.itemById("info").formattedText = _info_text(inputs, self.by_section)
        except Exception:
            _ui.messageBox("FoilBuilder:\n" + traceback.format_exc())


class BuildValidateHandler(adsk.core.ValidateInputsEventHandler):
    def __init__(self, by_section=False):
        super().__init__()
        self.by_section = by_section

    def notify(self, args):
        try:
            args.areInputsValid = _check(args.firingEvent.sender.commandInputs, self.by_section)[1] is None
        except Exception:
            args.areInputsValid = False


def _built_message(name, s):
    msg = f"Built {name}: root chord {s['c_root']:.1f} mm, AR {s['AR']:.1f}."
    if "strut_floor" in s:
        msg += (f"\nVertical: {s['strut_t']:.1f} mm thick, slot bottom at {s['strut_z_bot']:.1f} mm, "
                f"{s['strut_floor']:.1f} mm of foil under the slot.")
    if "rim_round_failed" in s:
        msg += ("\nThe slot rim round was skipped (Fusion could not make it); the rim is sharp. "
                f"Try a smaller rim round.\nFusion said: {s['rim_round_failed']}")
    return msg + ("\nEdit the foil_* parameters (Modify > Change Parameters) and use "
                  "'Rebuild Foil' to regenerate.")


class BuildExecuteHandler(adsk.core.CommandEventHandler):
    def __init__(self, by_section=False):
        super().__init__()
        self.by_section = by_section

    def notify(self, args):
        try:
            inputs = args.command.commandInputs
            p = _params_from_inputs(inputs)
            up = inputs.itemById("up_axis").selectedItem.name
            replace = inputs.itemById("replace").value
            foil_src, strut_src = _sources_from_inputs(inputs, self.by_section)
            name, s = build_foil(p, foil_src, up, replace, strut_src)
            _ui.messageBox(_built_message(name, s), APP_NAME)
        except Exception:
            _ui.messageBox("FoilBuilder build failed:\n" + traceback.format_exc())


class RebuildCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def notify(self, args):
        try:
            h = RebuildExecuteHandler(); args.command.execute.add(h); _handlers.append(h)
        except Exception:
            _ui.messageBox("FoilBuilder:\n" + traceback.format_exc())


class RebuildExecuteHandler(adsk.core.CommandEventHandler):
    def notify(self, args):
        try:
            design = _design()
            p, path, up = _read_params(design)
            secs = _read_sections(design)
            if secs:
                foil_src, strut_src = secs, _read_sections(design, "strut_sections")
            else:
                foil_src, strut_src = path, _read_attr(design, "strut_dat_path")
            files = [s for _, s in foil_src] if secs else [foil_src]
            if p["strut_on"]:
                files += [s for _, s in strut_src] if secs else [strut_src]
            missing = [f or "(none)" for f in files if not f or not os.path.isfile(f)]
            if missing:
                _ui.messageBox("Airfoil file(s) not found:\n" + "\n".join(missing) +
                               "\nRun 'Foil from .dat' or 'Foil from .dat by section' again.", APP_NAME)
                return
            build_foil(p, foil_src, up, True, strut_src)
        except Exception:
            _ui.messageBox("FoilBuilder rebuild failed:\n" + traceback.format_exc())


# ============================================================ add-in entry points
def _add_button(cmd_id, name, tooltip, handler):
    cdef = _ui.commandDefinitions.itemById(cmd_id)
    if cdef:
        cdef.deleteMe()
    cdef = _ui.commandDefinitions.addButtonDefinition(cmd_id, name, tooltip, "")
    cdef.commandCreated.add(handler)
    _handlers.append(handler)
    panel = _ui.workspaces.itemById(WORKSPACE_ID).toolbarPanels.itemById(PANEL_ID)
    ctrl = panel.controls.itemById(cmd_id)
    if ctrl:
        ctrl.deleteMe()
    panel.controls.addCommand(cdef)


def run(context):
    global _app, _ui
    try:
        _app = adsk.core.Application.get()
        _ui = _app.userInterface
        _add_button(CMD_BUILD, "Foil from .dat",
                    "Build a parametric T-foil lifting surface from an airfoil .dat file",
                    BuildCreatedHandler())
        _add_button(CMD_BUILD_SEC, "Foil from .dat by section",
                    "Build the foil from several airfoil .dat/.csv files placed along the span, "
                    "blending the section between them",
                    BuildCreatedHandler(by_section=True))
        _add_button(CMD_REBUILD, "Rebuild Foil",
                    "Rebuild the foil from the foil_* user parameters", RebuildCreatedHandler())
    except Exception:
        if _ui:
            _ui.messageBox("FoilBuilder failed to start:\n" + traceback.format_exc())


def stop(context):
    try:
        panel = _ui.workspaces.itemById(WORKSPACE_ID).toolbarPanels.itemById(PANEL_ID)
        for cid in (CMD_BUILD, CMD_BUILD_SEC, CMD_REBUILD):
            c = panel.controls.itemById(cid)
            if c:
                c.deleteMe()
            d = _ui.commandDefinitions.itemById(cid)
            if d:
                d.deleteMe()
    except Exception:
        if _ui:
            _ui.messageBox("FoilBuilder failed to stop:\n" + traceback.format_exc())
