---
title: FoilBuilder - Fusion add-in
tags: [aquafoil, fusion, cad, hydrofoil, tool]
created: 2026-10-19
---

# FoilBuilder (Fusion add-in) v0.2

Builds a hydrofoil lifting surface natively in Fusion from any airfoil `.dat` file:
- **Main foil**: bi-elliptical planform, with a flap split on a straight hinge line.
- **Rudder foil**: same planform, with curved ("banana") anhedral.

The add-in creates the construction planes, sketches, fitted splines, loft and split in Fusion, so the geometry is Fusion's own B-rep (no Python CAD libraries).

> [!warning] First version, not yet run inside Fusion
> The geometry maths (`foilgeom.py`) is tested: it reproduces the CadQuery Bonetti models exactly and reads Selig and Lednicer `.dat` files. The Fusion API calls have been checked for logic and plane/point consistency against a mock, but **not in a live Fusion session**. If something fails, the error dialog shows a traceback. Send it over and it can be fixed quickly.

## Install
1. Copy the whole `FoilBuilder` folder somewhere permanent. Keep the folder name, `FoilBuilder.py` and `FoilBuilder.manifest` together.
2. In Fusion: **Utilities → Add-Ins → Scripts and Add-Ins** (Shift+S) → **Add-Ins** tab → **+** (or the green plus), then select the `FoilBuilder` folder.
3. Select it and click **Run** (tick *Run on Startup* if you like).
4. Two buttons appear in **Solid → Create**: **Foil from .dat** and **Rebuild Foil**.

## Use
1. **Foil from .dat** opens a dialog:
   - **Browse...** to pick a `.dat` file (`airfoils/` has `e393.dat` and `NACA66-012.dat`). The info line shows the section's t/c and camber.
   - **Preset**: *Bonetti main foil* or *Bonetti rudder foil*, or set values by hand.
   - **Planform:** span (mm), area (cm²), Δ (ellipse axis position from the TE, fraction of root chord). The root chord is computed as 4A/(πb).
   - **Flap:** on/off, flap chord / root chord, minimum flap chord where the flap ends.
   - **Curved anhedral:** tip drop (0 = flat) and the ellipse vertical semi-axis.
   - **Advanced:** sections per half span, tip cut, TE thickness, vertical axis (Z or Y up), replace previous foil.
2. Click **OK**. The bodies are `foil`, or `main element` + `flap`, and their features are grouped in the timeline as **"FoilBuilder: Foil …"**.
   - In **Assembly** (or classic) documents the foil goes in its own component **"Foil …"**.
   - In **Part Design** documents (only one component allowed) it's built directly in the part.
3. All values are saved as **user parameters** (`foil_span`, `foil_area_cm2`, `foil_delta`, …). To change the design: **Modify → Change Parameters**, edit the values, then click **Rebuild Foil**. The previous foil component is replaced.

> [!note] Why "Rebuild" rather than live parameters
> A spline loaded from a `.dat` can't be tied to Fusion parameter expressions, so changing a parameter can't update the loft by itself. Rebuild Foil regenerates the geometry from the parameters in one click. Anything you model *on top* of the foil (bulb, strut, fillets) should go in a separate component, because it won't survive a rebuild.

## Coordinate system
- x aft (chordwise), y spanwise, z up (or Y up if selected, with span along -Z).
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
