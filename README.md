# stl2curves

Convert STL meshes (3D-printing downloads, exports from other programs) into solid
STEP files with **true curved surfaces**, so they can be edited properly in Fusion,
FreeCAD or any other CAD program.

A plain STL-to-STEP conversion gives you a solid made of thousands of tiny flat
triangles: you can't select a hole, measure a diameter or fillet an edge. stl2curves
works out which surface each part of the mesh really lies on and rebuilds it exactly:

| Rebuilt as | Examples |
|---|---|
| flat faces | walls, floors, straight chamfers (one face per flat area) |
| cylinders | holes, pins, bosses, straight rounded edges (fillets), slot ends |
| cones | countersinks, chamfers around round edges, cone tips |
| spheres | domes, dimples, ball ends, corners where three rounded edges meet |
| tori | rounded edges that follow a curve (rounded box corners, fillets round a pin) |

Surfaces can be cut to any outline (a hole running out through a sloped face, two
holes crossing), and overlapping bodies in one STL are fused automatically.

Areas that don't match any of these surfaces (freeform shapes, variable fillets)
stay as small flat facets, so the result always matches the mesh.

## Install

Python 3.10+:

```
pip install -r requirements.txt
```

## Use

```
python stl2curves.py part.stl                    # -> part.step next to it
python stl2curves.py some_folder --out STEP      # every .stl in a folder
python stl2curves.py a.stl b.stl --merge all     # also one combined STEP
python stl2curves.py part.stl --details          # list every rebuilt feature
python render.py part.step                       # PNG preview coloured by surface type
```

On Windows, `stl2curves.bat` accepts drag-and-drop (or put a shortcut to it in
`shell:sendto` for a right-click "Send to" entry).

Options: `--tol` sewing tolerance (mm), `--no-fuse` keep overlapping bodies separate,
`--no-curves` flat faces only.

## How it works

1. **Facets** — triangles are grouped into flat facets (coplanar pieces).
2. **Detection** (`features.py`) — neighbouring facets propose a surface (cylinder,
   cone, sphere; tori on axes already found), which is grown while mesh corners stay
   within a micron of it. Strict passes first (a patch must follow its surface's natural
   boundary lines, or meet its neighbours at creases), then a permissive pass for what
   is left, then neighbouring patches on the same surface are merged.
3. **Build** (`build.py`) — one exact face per patch and one planar face per flat
   area, with shared edges (exact lines/arcs where a patch's boundary follows its
   natural lines, splines through the mesh points otherwise), sewn into a solid.
4. **Checks** — every result must be a valid closed solid whose volume matches the
   mesh plus the predicted change from the curves; features that fail are left
   faceted instead of spoiling the part.

## Tests

```
python tests/run_tests.py
```

Builds five CAD test parts (holes, pins, chamfers, crossing holes, a filleted lid,
domes, balls...), meshes them, converts them back and checks each against the original
design (volume, validity, face count).

## Credits

The permissive final pass and the outline-polygon fallback face builder follow ideas
from [stlToSolid](https://github.com/Crypto69/stlToSolid) (no code copied). Built on
[OpenCascade](https://dev.opencascade.org/) via
[cadquery-ocp](https://pypi.org/project/cadquery-ocp/).
