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
2. In Fusion, open **Utilities → Add-Ins → Scripts and Add-Ins** (Shift+S), go to the **Add-Ins** tab, click **+** and select the `FoilBuilder` folder.
3. Select it and click **Run**. Tick *Run on Startup* if you want it every time.
4. Three buttons appear under **Solid → Create**: **Foil from .dat**, **Foil from .dat by section** and **Rebuild Foil**.

After updating the files, **Stop** and **Run** the add-in again so Fusion picks up the changes.

## Building a foil

**Foil from .dat** opens a dialog. The info line at the bottom summarises the result as you change values. If something can't be built, it gives the reason and **OK** stays disabled.

### Airfoil and preset
- **Browse...**: pick the airfoil (`.dat` or `.csv`, Selig or Lednicer layout). A CSV is `x,y` per row, and a header line is fine. `airfoils/` contains `e393.dat` and `NACA66-012.dat`. The info line shows the section's t/c and camber.
- **Preset**: *Bonetti main foil* or *Bonetti rudder foil* fill in the planform, flap, anhedral, bulb and advanced values. Choose *(keep current values)* to set them by hand. Presets don't touch the vertical support settings.

### Planform
- **Span** (mm) and **Area** (cm²). The root chord follows as 4A/(πb).
- **Δ axis from TE**: where the ellipse axis sits, as a fraction of the root chord measured from the trailing edge.

### Flap
- **Split flap** on/off, **Flap / root chord**, and **Flap min chord**: the flap ends where its chord drops below this.
- With a bulb, the flap is split into port and starboard halves round it.

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
- **Slot rim round** (default 0.5 mm, 0 = sharp): rounds the edge of the cutout, which gives the vertical a softer lead-in when you assemble it. With the blend on, this rounds the top of the collar instead (see below). If Fusion can't make the round, the foil is still built with a sharp rim and the message says so.

The info line shows the vertical's thickness, the deepest point of the slot, the least foil left under it, and any overhang. The vertical can't be built if it isn't over the foil at all, if its slot would overlap the flap, or if its top is too low.

### Blend into the vertical
Off by default. Tick **Blend foil up into the vertical** and the foil rises in a collar round the vertical before the cutout is made. The flow then meets the vertical closer to its angle instead of at a sharp 90° corner. The collar is part of the foil body, and the vertical still drops into the slot as a separate piece.

- **Blend height at the vertical** (default 3 mm): how far the collar stands above the foil where it meets the vertical. The rim follows the foil's top.
- **Blend width at the sides** (default 6 mm): how far out from the vertical's sides the collar reaches.
- **Blend length ahead of the vertical** and **behind the vertical** (default 6 mm each): how far the collar reaches forward of the vertical's nose and aft of its tail. Lengthen them for a longer, gentler ramp. They apply straight ahead and behind and taper smoothly into the side width, giving a teardrop shape. It gets a little wider than the sides just round the ends.
- **Blend angle at vertical** (default 25°): the collar's angle where it meets the vertical, measured from vertical. Smaller is closer to tangent. 0 would leave a knife edge at the rim.

The blend never blocks a build. It fades out smoothly where the foil gets thin (towards a nose or tail, or where the vertical overhangs) and before the flap.

With the blend on, **Slot rim round** rounds the collar's top edge, where it meets the vertical. That edge is a wedge as sharp as the blend angle, so only a small round fits. With the defaults (3 mm high, 25°) that's about 0.5 mm, and a larger blend angle allows more: about 1 mm at 45°. If the round is bigger than that, the info line warns you. If Fusion then can't make it, the build still finishes with a sharp edge.

How it's built: FoilBuilder works out the collar shape itself and lofts it like a fillet. It doesn't use Fusion's fillet tool, which can't handle this joint. The slices are cut at right angles to the vertical's wall all the way round it, fanning out round its nose and tail. Each has the same make-up and the same number of evenly spaced points, which keeps the surface smooth. Each slice reaches 1.2 mm under the foil's surface, and the collar meets the foil 8° off tangent. Both make the join a clean overlap rather than two surfaces just touching, which is what Fusion handles reliably. That hidden part is never visible, and where the collar fades out its surface sits 0.3 mm under the foil's.

### Advanced
- **Sections per half**: span stations on each side.
- **Tip cut**: the last station, as a fraction of the semi-span.
- **TE thickness**: trailing-edge thickness as % of local chord.
- **Vertical axis**: Y up (default) or Z up. See [Coordinate system](#coordinate-system).
- **Replace previous foil**: delete the previous FoilBuilder foil before building.

### The result
- **Bodies**: `foil`, or `main element` plus `flap` (`flap port` / `flap starboard` with a bulb). With the support on there is also a `vertical` body.
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
| Flap | `foil_flap_on`, `foil_flap_frac`, `foil_flap_min_chord` |
| Anhedral | `foil_anhedral_depth`, `foil_anhedral_semiaxis` |
| Advanced | `foil_n_half`, `foil_tip_eta`, `foil_te_thick_pct` |
| Bulb | `foil_bulb_on`, `foil_bulb_length`, `foil_bulb_nose_frac`, `foil_bulb_thickness`, `foil_bulb_half_width` |
| Vertical | `foil_strut_on`, `foil_strut_chord`, `foil_strut_le_x`, `foil_strut_height`, `foil_strut_slot_depth`, `foil_strut_rim_round` |
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
| `foilgeom.py` | Pure-Python geometry core (planform, sections, vertical, blend). Also usable outside Fusion, e.g. in optimisation scripts or checks |
| `FoilBuilder.manifest` | Add-in manifest |
| `airfoils/` | Example sections |

## Extending
New planform options (chord at N stations, twist, sweep) go in `foilgeom.Planform`, and the Fusion side doesn't need to change. Because the same `foilgeom.py` can run inside an optimisation loop, the optimiser and the CAD always describe the same geometry.
