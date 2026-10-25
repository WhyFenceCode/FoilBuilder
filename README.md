---
title: FoilBuilder - Fusion add-in
tags: [aquafoil, fusion, cad, hydrofoil, tool]
created: 2026-10-19
---

# FoilBuilder (Fusion add-in) v0.3

Builds a hydrofoil lifting surface natively in Fusion from any airfoil `.dat` file, or from several `.dat` / `.csv` airfoils placed along the span and blended between them:
- **Main foil**: bi-elliptical planform, with a flap split on a straight hinge line.
- **Rudder foil**: same planform, with curved ("banana") anhedral.
- Optional blended bulb at the centreline.
- Optional vertical support (strut): a straight, separate piece that slots into a cutout in the top of the foil/bulb.

The add-in creates the construction planes, sketches, fitted splines, loft and split in Fusion, so the geometry is Fusion's own B-rep (no Python CAD libraries).

> [!warning] First version, not yet run inside Fusion
> The geometry maths (`foilgeom.py`) is tested: it reproduces the CadQuery Bonetti models exactly and reads Selig and Lednicer `.dat` files. The Fusion API calls have been checked for logic and plane/point consistency against a mock, but **not in a live Fusion session**. If something fails, the error dialog shows a traceback. Send it over and it can be fixed quickly.

## Install
1. Copy the whole `FoilBuilder` folder somewhere permanent. Keep the folder name, `FoilBuilder.py` and `FoilBuilder.manifest` together.
2. In Fusion: **Utilities → Add-Ins → Scripts and Add-Ins** (Shift+S) → **Add-Ins** tab → **+** (or the green plus), then select the `FoilBuilder` folder.
3. Select it and click **Run** (tick *Run on Startup* if you like).
4. Three buttons appear in **Solid → Create**: **Foil from .dat**, **Foil from .dat by section** and **Rebuild Foil**.

## Use
1. **Foil from .dat** opens a dialog:
   - **Browse...** to pick a `.dat` file (`airfoils/` has `e393.dat` and `NACA66-012.dat`). The info line shows the section's t/c and camber.
   - **Preset**: *Bonetti main foil* or *Bonetti rudder foil*, or set values by hand.
   - **Planform:** span (mm), area (cm²), Δ (ellipse axis position from the TE, fraction of root chord). The root chord is computed as 4A/(πb).
   - **Flap:** on/off, flap chord / root chord, minimum flap chord where the flap ends.
   - **Curved anhedral:** tip drop (0 = flat) and the ellipse vertical semi-axis.
   - **Bulb (blended):** on/off, length at the centreline, share of the extra length ahead of the LE, max thickness, and the spanwise half-width over which it blends into the wing. With a flap, the flap is split into port and starboard halves around the bulb.
   - **Vertical support:** off by default. Tick it to add a vertical (strut) on the centreline, then pick its airfoil `.dat` / `.csv` (chord along x, thickness spanwise) and set:
     - **Chord** and **LE aft of foil nose**: chordwise position, measured from the foil's nose on the centreline (the bulb nose when there is a bulb).
     - **Top height above chord line**: the vertical is straight, with the same section all the way from the slot bottom to a flat top face at this height, measured from the foil chord line at the centreline. Extrude from that face to lengthen it.
     - **Slot depth into foil**: the vertical reaches this far below the lowest point of the foil top under it, ending in a flat bottom. The same shape is cut out of the foil as the slot, so the vertical drops in and can be screwed from below. Screw holes aren't modelled.
     - **Slot rim round**: a small round on the edge of the cutout where the vertical enters the foil (default 1 mm, 0 = sharp). If Fusion can't make it, the foil is still built with a sharp rim, and the message says so.
     - The vertical's section is closed with a tiny rounded TE (at least 0.5 mm thick), so the slot rim has no sharp corners.
     - The info line shows the vertical's thickness and how much foil is left under the slot (at least 2 mm). If it doesn't fit (slot too deep, vertical outside the foil chord, slot overlapping the flap, or top below the foil), the reason is shown and **OK** stays disabled.
   - **Advanced:** sections per half span, tip cut, TE thickness, vertical axis (Y up by default, or Z up), replace previous foil.
2. Click **OK**. The bodies are `foil`, or `main element` + `flap`, plus a separate `vertical` body when the support is on, and their features are grouped in the timeline as **"FoilBuilder: Foil …"**.
   - In **Assembly** (or classic) documents the foil goes in its own component **"Foil …"**.
   - In **Part Design** documents (only one component allowed) it's built directly in the part.
3. **Foil from .dat by section** is the same dialog, but the section can change along the span:
   - **Add...** picks one or more `.dat` / `.csv` files (Selig or Lednicer layout; a CSV is `x,y` per row, a header line is fine). Each file becomes a row with its position **η** = fraction of the semi-span (0 = centreline, 1 = tip), mirrored to both sides. Picking several files into an empty table spreads them evenly from root to tip; edit η in the table.
   - **Remove** deletes the selected row (or the last one).
   - Every airfoil is resampled onto the same cosine-spaced x stations, and each span station gets a linear blend of the two neighbouring rows. Inboard of the first row / outboard of the last, that row's airfoil is used unchanged. A section is placed exactly at each row's η.
   - The **Vertical support** group gets its own section table. Here η is the height fraction, from the highest point of the foil top under the vertical (0) to its flat top (1). The tenon uses the η = 0 section.
   - The file list and positions are stored in the design, so **Rebuild Foil** rebuilds by section too. Running plain **Foil from .dat** switches the design back to a single section.
4. All values are saved as **user parameters** (`foil_span`, `foil_area_cm2`, `foil_delta`, …). To change the design: **Modify → Change Parameters**, edit the values, then click **Rebuild Foil**. The previous foil component is replaced.

> [!note] Why "Rebuild" rather than live parameters
> A spline loaded from a `.dat` can't be tied to Fusion parameter expressions, so changing a parameter can't update the loft by itself. Rebuild Foil regenerates the geometry from the parameters in one click. Anything you model *on top* of the foil (bulb, strut, fillets) should go in a separate component, because it won't survive a rebuild.

## Coordinate system
- Default (Y up): model X aft (chordwise), Y up, span along -Z. With Z up selected: X aft, Y spanwise, Z up.
- The choice is stored in the design, so Rebuild Foil keeps it; designs built before with Z up stay Z up until you pick Y in the dialog.
- Origin: on the planform's ellipse axis at the centreline.
- Foil incidence is 0°. Set it at assembly (strut incidence / rudder rake).

## Files
| File | Purpose |
|---|---|
| `FoilBuilder.py` | Fusion add-in (UI, Fusion geometry calls) |
| `foilgeom.py` | Pure-Python geometry core. Also usable outside Fusion (optimization scripts, checks) |
| `FoilBuilder.manifest` | Add-in manifest |
| `airfoils/` | Example sections |

## Extending
New planform options go in `foilgeom.Planform` (e.g. chord at N stations, twist, sweep), and the Fusion side doesn't change. Because the same `foilgeom.py` can run in the optimization loop, the optimizer and the CAD always describe the same geometry.
