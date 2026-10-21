"""FoilBuilder: Fusion add-in that builds a parametric T-foil (main or rudder
lifting foil) from any airfoil .dat file.

Buttons (Solid > Create panel):
  * Foil from .dat   - dialog: pick .dat, preset, parameters -> builds the foil
  * Rebuild Foil     - reads the foil_* user parameters and rebuilds

Geometry maths lives in foilgeom.py (pure Python). This file only talks to
Fusion: construction planes, sketches, fitted splines, loft, flap split.
Fusion's internal length unit is cm; foilgeom works in mm (MM2CM below).
"""
import adsk.core
import adsk.fusion
import os
import sys
import math
import time
import traceback

_HERE = os.path.dirname(os.path.realpath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import foilgeom as G  # noqa: E402

APP_NAME = "FoilBuilder"
CMD_BUILD = "aquafoil_FoilBuilder_build"
CMD_REBUILD = "aquafoil_FoilBuilder_rebuild"
PANEL_ID = "SolidCreatePanel"
WORKSPACE_ID = "FusionSolidEnvironment"
ATTR_GROUP = "aQuaFoil_FoilBuilder"
MM2CM = 0.1

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
]
BOOL_KEYS = ("flap_on", "bulb_on")

_app = None
_ui = None
_handlers = []


# ============================================================ helpers
def _design():
    d = adsk.fusion.Design.cast(_app.activeProduct)
    if not d:
        raise RuntimeError("Open a Fusion design (Design workspace) first.")
    return d


def _write_params(design, p, dat_path, up_axis):
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
    design.attributes.add(ATTR_GROUP, "dat_path", dat_path)
    design.attributes.add(ATTR_GROUP, "up_axis", up_axis)


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
    return p, (a.value if a else ""), (u.value if u else "Z")


def _default_up_axis():
    try:
        pref = _app.preferences.generalPreferences.defaultModelingOrientation
        if pref == adsk.core.DefaultModelingOrientations.YUpModelingOrientation:
            return "Y"
        if pref == adsk.core.DefaultModelingOrientations.ZUpModelingOrientation:
            return "Z"
        if pref == adsk.core.DefaultModelingOrientations.XUpModelingOrientation:
            return "X"
    except Exception:
        pass
    return "Y"


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
def build_foil(p, dat_path, up_axis="Z", replace=True):
    design = _design()
    with open(dat_path, "r", errors="ignore") as f:
        af_name, af = G.parse_dat(f.read())
    pf = G.Planform(p)
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
        result = _build_geometry(design, comp, af, pf, P)
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

    _write_params(design, p, dat_path, up_axis)
    return af_name, pf.summary()


def _build_geometry(design, comp, af, pf, P):
    helper = _tag(comp.sketches.add(comp.xYConstructionPlane))
    helper.name = "FoilBuilder helper points"

    profiles = []
    for y in pf.stations():
        z0, dz = pf.z_span(y)
        phi = math.atan(dz)
        ny, nz = -math.sin(phi), math.cos(phi)
        plane = _plane_through(comp, helper, P(0.0, y, z0), P(10.0, y, z0),
                               P(0.0, y + 10.0 * ny, z0 + 10.0 * nz))
        sk = _tag(comp.sketches.add(plane))
        sk.name = f"section y={y:.1f}"
        sk.isComputeDeferred = True
        coll = adsk.core.ObjectCollection.create()
        for (x, yy, zz) in pf.section(af, y):
            coll.add(sk.modelToSketchSpace(P(x, yy, zz)))
        spl = sk.sketchCurves.sketchFittedSplines.add(coll)
        sk.sketchCurves.sketchLines.addByTwoPoints(spl.endSketchPoint, spl.startSketchPoint)
        sk.isComputeDeferred = False
        if sk.profiles.count != 1:
            raise RuntimeError(f"Section at y={y:.1f} mm gave {sk.profiles.count} profiles")
        profiles.append(sk.profiles.item(0))
        sk.isVisible = False
        adsk.doEvents()
    helper.isVisible = False

    lofts = comp.features.loftFeatures
    li = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    for pr in profiles:
        li.loftSections.add(pr)
    li.isSolid = True
    loft = _tag(lofts.add(li))
    foil = loft.bodies.item(0)
    foil.name = "foil"

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
    return True


# ============================================================ dialog
def _info_text(inputs):
    try:
        p = _params_from_inputs(inputs)
        path = inputs.itemById("dat_path").text
        lines = []
        if path and os.path.isfile(path):
            with open(path, "r", errors="ignore") as f:
                name, af = G.parse_dat(f.read())
            t, xt, cm, xc = G.airfoil_stats(af)
            lines.append(f"<b>{name}</b>: {len(af)} pts, t/c {t*100:.1f}% @ {xt*100:.0f}%, "
                         f"camber {cm*100:.1f}% @ {xc*100:.0f}%")
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
        return "<br>".join(lines)
    except Exception as e:
        return f"<i>{e}</i>"


def _params_from_inputs(inputs):
    g = lambda i: inputs.itemById(i)
    return dict(
        span=g("span").value, area_cm2=g("area_cm2").value, delta=g("delta").value,
        flap_on=g("flap_on").value, flap_frac=g("flap_frac").value,
        flap_min_chord=g("flap_min_chord").value, anhedral_depth=g("anhedral_depth").value,
        anhedral_semiaxis=g("anhedral_semiaxis").value, n_half=g("n_half").value,
        tip_eta=g("tip_eta").value, te_thick_pct=g("te_thick_pct").value,
        bulb_on=g("bulb_on").value, bulb_length=g("bulb_length").value,
        bulb_nose_frac=g("bulb_nose_frac").value, bulb_thickness=g("bulb_thickness").value,
        bulb_half_width=g("bulb_half_width").value)


def _set_inputs(inputs, p):
    for key in ("span", "area_cm2", "delta", "flap_frac", "flap_min_chord", "anhedral_depth",
                "anhedral_semiaxis", "tip_eta", "te_thick_pct", "bulb_length", "bulb_nose_frac",
                "bulb_thickness", "bulb_half_width"):
        inputs.itemById(key).value = float(p.get(key, G.DEFAULTS[key]))
    inputs.itemById("n_half").value = int(p["n_half"])
    inputs.itemById("flap_on").value = bool(p["flap_on"])
    inputs.itemById("bulb_on").value = bool(p.get("bulb_on", False))


class BuildCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def notify(self, args):
        try:
            cmd = args.command
            inputs = cmd.commandInputs
            try:
                p, dat, up = _read_params(_design())
            except Exception:
                p, dat, up = dict(G.DEFAULTS), "", _default_up_axis()

            btn = inputs.addBoolValueInput("browse", "Airfoil .dat", False, "", False)
            try:
                btn.text = "Browse..."
            except Exception:
                pass
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

            g = inputs.addGroupCommandInput("g_adv", "Advanced").children
            g.addIntegerSpinnerCommandInput("n_half", "Sections per half", 4, 60, 1, int(p["n_half"]))
            g.addFloatSpinnerCommandInput("tip_eta", "Tip cut (fraction of semi-span)", "", 0.9, 0.9999, 0.001, p["tip_eta"])
            g.addFloatSpinnerCommandInput("te_thick_pct", "TE thickness (% chord)", "", 0.05, 3, 0.05, p["te_thick_pct"])
            ud = g.addDropDownCommandInput("up_axis", "Vertical axis", adsk.core.DropDownStyles.TextListDropDownStyle)
            ud.listItems.add("Z", up != "Y", "")
            ud.listItems.add("Y", up == "Y", "")
            g.addBoolValueInput("replace", "Replace previous foil", True, "", True)

            inputs.addTextBoxCommandInput("info", "", _info_text(inputs), 3, True)

            h1 = BuildExecuteHandler(); cmd.execute.add(h1); _handlers.append(h1)
            h2 = BuildChangedHandler(); cmd.inputChanged.add(h2); _handlers.append(h2)
            h3 = BuildValidateHandler(); cmd.validateInputs.add(h3); _handlers.append(h3)
        except Exception:
            _ui.messageBox("FoilBuilder dialog failed:\n" + traceback.format_exc())


class BuildChangedHandler(adsk.core.InputChangedEventHandler):
    def notify(self, args):
        try:
            inputs = args.firingEvent.sender.commandInputs   # top level (args.inputs is per group)
            changed = args.input
            if changed.id == "browse":
                dlg = _ui.createFileDialog()
                dlg.title = "Select airfoil coordinate file"
                dlg.filter = "Airfoil (*.dat;*.txt);;All files (*.*)"
                if dlg.showOpen() == adsk.core.DialogResults.DialogOK:
                    inputs.itemById("dat_path").text = dlg.filename
            elif changed.id == "preset":
                name = changed.selectedItem.name if changed.selectedItem else ""
                if name in G.PRESETS:
                    _set_inputs(inputs, G.PRESETS[name])
            inputs.itemById("info").formattedText = _info_text(inputs)
        except Exception:
            _ui.messageBox("FoilBuilder:\n" + traceback.format_exc())


class BuildValidateHandler(adsk.core.ValidateInputsEventHandler):
    def notify(self, args):
        try:
            inputs = args.firingEvent.sender.commandInputs
            path = inputs.itemById("dat_path").text
            ok = bool(path) and os.path.isfile(path)
            if ok:
                G.Planform(_params_from_inputs(inputs))
            args.areInputsValid = ok
        except Exception:
            args.areInputsValid = False


class BuildExecuteHandler(adsk.core.CommandEventHandler):
    def notify(self, args):
        try:
            inputs = args.command.commandInputs
            p = _params_from_inputs(inputs)
            path = inputs.itemById("dat_path").text
            up = inputs.itemById("up_axis").selectedItem.name
            replace = inputs.itemById("replace").value
            name, s = build_foil(p, path, up, replace)
            _ui.messageBox(f"Built {name}: root chord {s['c_root']:.1f} mm, AR {s['AR']:.1f}.\n"
                           "Edit the foil_* parameters (Modify > Change Parameters) and use "
                           "'Rebuild Foil' to regenerate.", APP_NAME)
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
            if not path or not os.path.isfile(path):
                _ui.messageBox("No valid .dat path stored in this design. Run 'Foil from .dat' first.", APP_NAME)
                return
            name, s = build_foil(p, path, up, True)
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
        _add_button(CMD_REBUILD, "Rebuild Foil",
                    "Rebuild the foil from the foil_* user parameters", RebuildCreatedHandler())
    except Exception:
        if _ui:
            _ui.messageBox("FoilBuilder failed to start:\n" + traceback.format_exc())


def stop(context):
    try:
        panel = _ui.workspaces.itemById(WORKSPACE_ID).toolbarPanels.itemById(PANEL_ID)
        for cid in (CMD_BUILD, CMD_REBUILD):
            c = panel.controls.itemById(cid)
            if c:
                c.deleteMe()
            d = _ui.commandDefinitions.itemById(cid)
            if d:
                d.deleteMe()
    except Exception:
        if _ui:
            _ui.messageBox("FoilBuilder failed to stop:\n" + traceback.format_exc())
