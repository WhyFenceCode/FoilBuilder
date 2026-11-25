---
title: FoilBuilder - Fusion add-in
tags: [aquafoil, fusion, cad, hydrofoil, tool]
created: 2026-10-19
---

# FoilBuilder (Fusion add-in) v0.4

FoilBuilder builds a hydrofoil lifting surface natively in Fusion from airfoil coordinate files. It creates the construction planes, sketches, fitted splines, lofts and splits itself, so the result is ordinary Fusion geometry with no external CAD libraries.

What it can build:
- **Main foil**: bi-elliptical planform, optionally with a flap split along a straight hinge line.
- **Rudder foil**: the same planform with curved ("banana") anhedral.
- **Varying section** (optional): several airfoils placed along the span, with the section blended between them.
- **Blended bulb** (optional): a bulb at the centreline.
- **Vertical support** (optional): a straight, separate strut that slots into a cutout in the top of the foil or bulb. The foil can rise in a smooth collar round it.

> [!note] Status
> The geometry maths in `foilgeom.py` is tested on its own: it reproduces the CadQuery Bonetti models and reads Selig and Lednicer files. The Fusion side has been run in Fusion during development, but the newest parts (the vertical's blend in particular) have had little use there. If a build fails, the error dialog shows a traceback. Send it over and it can be fixed quickly.

## Install

1. Copy the whole `FoilBuilder` folder somewhere permanent. Keep `FoilBuilder.py`, `foilgeom.py` and `FoilBuilder.manifest` together in a folder named `FoilBuilder`.
2. Install numpy into the add-in's own `lib` folder (see [numpy](#numpy) below). Fusion's built-in Python doesn't include it.
3. In Fusion, open **Utilities → Add-Ins → Scripts and Add-Ins** (Shift+S), go to the **Add-Ins** tab, click **+** and select the `FoilBuilder` folder.
4. Select it and click **Run**. Tick *Run on Startup* if you want it every time.
5. Three buttons appear under **Solid → Create**: **Foil from .dat**, **Foil from .dat by section** and **Rebuild Foil**.

After updating the files, **Stop** and **Run** the add-in again so Fusion picks up the changes.

### numpy
The geometry code uses numpy. Fusion's Python ships without numpy and without pip, so numpy goes into a `lib` folder inside `FoilBuilder`, which the add-in adds to its path. That numpy has to be built for the same Python version as Fusion's. At the time of writing that's Python 3.14, 64-bit; check `python3XX.dll` in Fusion's install folder under `%LOCALAPPDATA%\Autodesk\webdeploy\production\`.

With a standalone Python of that version installed, run this from the `FoilBuilder` folder:

```bash
py -3.14 -m pip install --target lib numpy
```

`lib/` is about 55 MB and is kept out of git (`.gitignore`). If numpy is missing or doesn't match Fusion's Python, the add-in shows a message saying so when it starts, instead of adding its buttons. When an update moves Fusion to a newer Python, delete `lib` and run the command again with the new version.

## Building a foil

**Foil from .dat** opens a dialog. The info line at the bottom summarises the result as you change values. If something can't be built, it gives the reason and **OK** stays disabled.

### Airfoil and preset
- **Browse...**: pick the airfoil (`.dat` or `.csv`, Selig or Lednicer layout). A CSV is `x,y` per row, and a header line is fine. `airfoils/` contains `e393.dat` and `NACA66-012.dat`. The info line shows the section's t/c and camber.
- **Preset**: *Bonetti main foil* or *Bonetti rudder foil* fill in the planform, flap, anhedral, bulb and advanced values. Choose *(keep current values)* to set them by hand. Presets don't touch the vertical support settings.

### Planform
- **Span** (mm) and **Area** (cm²). The root chord follows as 4A/(πb).
- **Δ axis from TE**: where the ellipse axis sits, as a fraction of the root chord measured from the trailing edge.

### Flap
- **Split flap** on/off and **Flap / root chord**: where the straight hinge line sits, as the flap's chord at the centreline.
- The flap is one piece behind the straight hinge line. With a bulb, the hinge line runs straight through it, so the flap carries the bulb's tail. Outboard it runs out until it tapers to nothing, where the trailing edge curves forward across the hinge line (about 95% of the semi-span with the defaults). The info line gives the flap's span. Its outer ends get very thin, down to the trailing-edge thickness.
- **Flap angle** (default 0°, ±45°, + = trailing edge down): the flap is built deflected by this. It turns as a rigid piece about a straight hinge axis along the hinge line, at the top of the joint between main element and flap on the centreline. Away from the centreline the joint's top edge sits a little lower (the wing is thinner than the bulb's crown), so out there the flap's front edge moves slightly fore or aft as it turns. Turning down, the flap's front face swings forward below the pivot into the main element (the two bodies overlap there); turning up, a wedge opens on the underside.
- **Join flap into the foil body**: after the split (and the turn), the flap is joined back into the main element as one `foil` body. A deflected flap only touches the main element along the hinge, so if Fusion can't join it, the flap is left as its own body and the message says so.

### Curved anhedral
- **Tip drop** (0 = flat) and the **Ellipse semi-axis**, the vertical semi-axis of the ellipse the span line follows.

### Bulb (blended)
- On/off, **Length at centreline**, **Share of extra length ahead of LE**, **Max thickness**, and **Blend half-width**: how far either side of the centreline the bulb blends into the wing.

### Vertical support
Off by default. Tick **Vertical support** to add a vertical (strut) on the centreline. It is its own straight body, with the same section from bottom to top, and slots into a matching cutout in the foil. It can then be screwed in from below; screw holes aren't modelled.

- **Browse...**: the vertical's airfoil, with chord along x and thickness across the span. Its trailing edge is rounded off at no less than 0.5 mm thick, so the cutout has no sharp corners.
- **Chord** (default 80 mm).
- **LE aft of foil nose** (default 25 mm): measured on the centreline from the foil's nose, which is the bulb nose when there's a bulb. Negative values put the LE ahead of the nose. The vertical can sit anywhere along the foil, even overhanging its nose or tail. Where it overhangs, its bottom is trimmed level with that end of the foil, so nothing hangs below.
- **Top height above chord line** (default 150 mm): height of the vertical's flat top face, measured from the foil's chord line at the centreline. Extrude that face to lengthen the vertical.
- **Slot depth into foil** (default 5 mm): the vertical's bottom, and the cutout, follow the foil's top surface this far down, always leaving at least 2 mm of foil underneath. Where the foil is too thin for that, the vertical simply sits on the surface. If the vertical overhangs a blunt bulb nose or tail, the last few mm can leave less than 2 mm, because the trim drops to that end's height so there's no gap at the vertical's sides. The info line reports the actual minimum.
  - **0** means no slot: the vertical sits on the foil, and its bottom follows the foil's top surface exactly, including the crown of the bulb across its thickness. The foil isn't cut. With the blend on, the collar is still cut back where the vertical stands, down to the foil's own top, so the vertical sits flush inside the collar.
- **Slot rim round** (default 0.5 mm, 0 = sharp): rounds the edge of the cutout, which gives the vertical a softer lead-in when you assemble it. With the blend on, this rounds the top of the collar instead (see below). With a slot depth of 0 and no blend there's no cutout, so nothing is rounded. If Fusion can't make the round, the foil is still built with a sharp rim and the message says so.

- **Join vertical into the foil body**: the vertical is joined into the foil (the main element when there's a flap) at the end, instead of staying a separate piece. No slot is cut and no rim round is made, since the two become one solid.
- **Blend foil up into the vertical**: see below. Leave it unticked for a plain vertical in a plain slot.

#### Flap relief
When the vertical reaches back over the flap, a cutout is made in its back so the flap can deflect up into that space. FoilBuilder takes the flap's profile under the vertical, rotates it about the hinge from 0° up to the maximum angle, and cuts away everything under that swept shape plus the clearance. The slot in the foil then stops ahead of the hinge, so the vertical only holds on the main element and the flap stays free. If the vertical ends ahead of the flap, nothing is cut.

- **Flap relief: max flap up** (default 15°): the largest upward deflection the vertical has to clear. The flap pivots about the top of the joint on the centreline (see [Flap](#flap)). If the flap is built set further up than this, the cutout clears that angle instead.
- **Flap relief: clearance** (default 1 mm): gap kept between the swept flap and the vertical.
- **Flap relief: front round** (default 5 mm, 0 = sharp): the front of the cutout is a quarter circle of this radius, reaching forward from the hinge, so there's no sharp inside corner in the vertical.

The info line shows the vertical's thickness, the deepest point of the slot, the least foil left under it, any overhang, and the size of the flap relief. The vertical can't be built if it isn't over the foil at all, if the flap relief would cut away its front, or if its top is too low.

### Blend into the vertical
Off by default. Tick **Blend foil up into the vertical** and the foil rises in a collar round the vertical before the cutout is made. The flow then meets the vertical closer to its angle instead of at a sharp 90° corner. The collar is part of the foil body, and the vertical still drops into the slot as a separate piece.

- **Blend height at the vertical** (default 3 mm): how far the collar stands above the foil where it meets the vertical. The rim follows the foil's top.
- **Blend width at the sides** (default 6 mm): how far out from the vertical's sides the collar reaches.
- **Blend length ahead of the vertical** and **behind the vertical** (default 6 mm each): how far the collar reaches forward of the vertical's nose and aft of its tail. Lengthen them for a longer, gentler ramp. They apply straight ahead and behind and taper smoothly into the side width, giving a teardrop shape. It gets a little wider than the sides just round the ends.
- **Blend angle at vertical** (default 25°): the collar's angle where it meets the vertical, measured from vertical. Smaller is closer to tangent. 0 would leave a knife edge at the rim.

The blend never blocks a build. It fades out smoothly where the foil gets thin (towards a nose or tail, or where the vertical overhangs) and before the flap (or before the vertical's flap relief).

With the blend on, **Slot rim round** rounds the collar's top edge, where it meets the vertical. That edge is a wedge as sharp as the blend angle, so only a small round fits. With the defaults (3 mm high, 25°) that's about 0.5 mm, and a larger blend angle allows more: about 1 mm at 45°. If the round is bigger than that, the info line warns you. If Fusion then can't make it, the build still finishes with a sharp edge.

How it's built: FoilBuilder works out the collar shape itself and lofts it like a fillet. It doesn't use Fusion's fillet tool, which can't handle this joint. The slices are cut at right angles to the vertical's wall all the way round it, fanning out round its nose and tail. Each has the same make-up and the same number of evenly spaced points, which keeps the surface smooth. Each slice reaches 1.2 mm under the foil's surface, and the collar meets the foil 8° off tangent. Both make the join a clean overlap rather than two surfaces just touching, which is what Fusion handles reliably. That hidden part is never visible, and where the collar fades out its surface sits 0.3 mm under the foil's.

### Advanced
- **Sections per half**: span stations on each side.
- **Tip cut**: the last station, as a fraction of the semi-span.
- **TE thickness**: trailing-edge thickness as % of local chord.
- **Vertical axis**: Y up (default) or Z up. See [Coordinate system](#coordinate-system).
- **Replace previous foil**: delete the previous FoilBuilder foil before building.

### The result
- **Bodies**: `foil`, or `main element` plus `flap`. With the support on there is also a `vertical` body, unless it's joined into the foil. With the flap joined there's a single `foil` body.
- **Timeline**: the features are grouped as **"FoilBuilder: Foil …"**.
- **Assembly** (or classic) documents: the foil goes in its own component, **"Foil …"**.
- **Part Design** documents (one component only): it's built directly in the part.

## Building by section

**Foil from .dat by section** is the same dialog, but the airfoil can change along the span. In place of the single file picker there is a table with one row per airfoil.

- **Add...**: pick one or more `.dat` / `.csv` files. Each becomes a row with its position **η**, the fraction of the semi-span (0 = centreline, 1 = tip), mirrored to both sides. Several files added to an empty table are spread evenly from root to tip; edit η in the table.
- **Remove selected (or last)**: deletes the selected row, or the last row if none is selected.
- **Blending**: every airfoil is resampled onto the same chordwise stations, and each span station gets a linear blend of its two neighbouring rows. Inboard of the first row and outboard of the last, that row's airfoil is used unchanged. A station is placed exactly at each row's η.
- **Vertical sections**: the Vertical support group gets its own table. There η is height, from the highest point of the foil top under the vertical (0) to the vertical's flat top (1). Everything below η = 0 uses the η = 0 section.

## Changing a design: the parameter table

Every dialog value is saved in the design as a Fusion **user parameter** with a `foil_` prefix. Values you can change there:

| Group | Parameters |
|---|---|
| Planform | `foil_span`, `foil_area_cm2`, `foil_delta` |
| Flap | `foil_flap_on`, `foil_flap_frac`, `foil_flap_angle`, `foil_join_flap` |
| Anhedral | `foil_anhedral_depth`, `foil_anhedral_semiaxis` |
| Advanced | `foil_n_half`, `foil_tip_eta`, `foil_te_thick_pct` |
| Bulb | `foil_bulb_on`, `foil_bulb_length`, `foil_bulb_nose_frac`, `foil_bulb_thickness`, `foil_bulb_half_width` |
| Vertical | `foil_strut_on`, `foil_strut_chord`, `foil_strut_le_x`, `foil_strut_height`, `foil_strut_slot_depth`, `foil_strut_rim_round`, `foil_join_strut` |
| Flap relief | `foil_strut_flap_up_deg`, `foil_strut_relief_clear`, `foil_strut_relief_round` |
| Blend | `foil_blend_on`, `foil_blend_height`, `foil_blend_width`, `foil_blend_front_width`, `foil_blend_back_width`, `foil_blend_wall_angle` |

On/off values are `1` or `0`. The airfoil files and their η positions, and the vertical axis, are stored in the design too, but not as parameters. To change those, run the dialog again.

To edit:
1. Open the table: **Solid → Modify → Change Parameters** (the *fx* icon). The `foil_` values are listed under **User Parameters**.
2. Change values in the **Expression** column and click **OK**.
3. Click **Rebuild Foil** (Solid → Create). It replaces the previous foil with one built from the new values.

> [!note] Why "Rebuild" rather than live parameters
> A spline loaded from a coordinate file can't be tied to Fusion parameter expressions, so changing a parameter can't update the loft by itself. Rebuild Foil regenerates everything in one click. Anything you model yourself on top of the foil (fittings, extra fillets, an extension of the vertical) should go in a separate component, because a rebuild doesn't keep it.

## Coordinate system
- **Y up** (default): model X points aft (chordwise), Y is up, and the span runs along −Z.
- **Z up**: X aft, Y spanwise, Z up.
- The choice is stored in the design, so Rebuild Foil keeps it. A design built with Z up stays Z up until you pick Y in the dialog.
- The origin is on the planform's ellipse axis at the centreline, on the chord line.
- The foil is built at 0° incidence. Set incidence at assembly (strut incidence or rudder rake).

## Files
| File | Purpose |
|---|---|
| `FoilBuilder.py` | Fusion add-in: dialogs and Fusion geometry calls |
| `foilgeom.py` | Geometry core in Python and numpy (planform, sections, vertical, blend). Also usable outside Fusion, e.g. in optimisation scripts or checks |
| `FoilBuilder.manifest` | Add-in manifest |
| `airfoils/` | Example sections |
| `lib/` | numpy for Fusion's Python (installed locally, not in git; see [numpy](#numpy)) |

## Extending
New planform options (chord at N stations, twist, sweep) go in `foilgeom.Planform`, and the Fusion side doesn't need to change. Because the same `foilgeom.py` can run inside an optimisation loop, the optimiser and the CAD always describe the same geometry.
