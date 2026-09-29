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
| smooth blends | where rounded edges meet at a corner, fillets along spline-shaped edges |

Surfaces can be cut to any outline (a hole running out through a sloped face, two
holes crossing), and overlapping bodies in one STL are fused automatically, as are
bodies that only touch along an edge.

Curved areas that none of the exact surfaces fit (corner blends, rounded edges
following a spline) become one smooth freeform face each, fitted through the mesh
corners and rolling tangentially into the faces beside them. Anything else stays as
flat facets, so the result always matches the mesh.

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
python stl2curves.py part.stl --true-size        # rebuild at the apparent design size
python render.py part.step                       # PNG preview coloured by surface type
```

On Windows, `stl2curves.bat` accepts drag-and-drop (or put a shortcut to it in
`shell:sendto` for a right-click "Send to" entry).

Options: `--tol` sewing tolerance (mm), `--no-fuse` keep overlapping bodies separate,
`--no-curves` flat faces only, `--no-blends` no smooth freeform faces (exact surfaces
and flat facets only).

### Design size

People design with round numbers (a 5 mm hole, a 0.5 mm fillet, a 1/4" slot), but an
STL is often off by one overall factor: scaled to 99% in a slicer for fit, exported in
inches or centimetres, and so on. Every conversion reports the size the part appears to
have been designed at: the single scale factor that turns the most of its radii and
wall thicknesses into round numbers (0.5 mm or 0.1 mm steps, 1/16" for inch designs),
e.g. `looks scaled to 99% (x1.01 gives round mm: 0.495 -> 0.5, 2.97 -> 3 ...)`.
`--true-size` rebuilds the part at that size with radii snapped to the round values.
Without it the size is left alone (a slight scale is often deliberate, for fit).

## How it works

1. **Facets** — triangles are grouped into flat facets (coplanar pieces). The mesh is
   tidied first: sliver triangles (three corners in a line) are split away, and when the
   file rounded its coordinates (many exports write 0.001 mm steps) the rounding noise is
   allowed for, both when deciding which triangles are coplanar (long thin triangles
   tilt noticeably) and in every surface test. Pieces of the mesh that merely touch
   (stacked blocks) are converted separately and fused.
2. **Detection** (`features.py`) — first, each smooth area bounded by sharp edges is
   tried as a whole (a countersink, a plain hole): one cylinder, cone or sphere fitted
   through all its corners at once, allowing a few outline corners that a neighbouring
   face's triangulation put slightly off the curve. Then neighbouring facets propose a
   surface (cylinder, cone, sphere; tori on axes already found), which is grown while
   mesh corners stay within a micron of it. Strict passes first (a patch must follow its surface's natural
   boundary lines, or meet its neighbours at creases), then fillets between flat faces
   (`fillets.py`: two faces meeting at an angle, rounded off, give a cylinder touching
   both, so only the radius is unknown and every strip corner pins it down; this finds
   fillets however coarsely or irregularly they were cut into triangles), then a
   permissive pass for what is left, then neighbouring patches on the same surface are
   merged, chains of short
   cylinder strips or sphere bits are replaced by the torus they approximate, and equal
   radii and near-axis-aligned axes are snapped to exact values (each snap kept only if
   the patch still lies on the mesh).
3. **Blends** (`blends.py`) — dense curved areas left over (unrecognised facets, and
   small pieces that only stand in for a surface they don't really fit) are grouped
   into regions that each turn no more than one smooth face can follow.
4. **Build** (`build.py`) — one exact face per patch, one smooth N-sided patch per
   blend (through the mesh corners; blends are built on their own first, a blend that
   won't fit is cut in two and retried) and one planar face per flat area, with shared
   edges (exact lines/arcs where a patch's boundary
   follows its natural lines, splines through the mesh points otherwise, split at sharp
   corners), sewn into a solid. Bodies that touch along an edge (`bodies.py`) are built
   separately and fused.
5. **Checks** — every result must be a valid closed solid whose volume matches the
   mesh plus the predicted change from the curves; features that fail are left
   faceted instead of spoiling the part (a blend that fails gives back the pieces it
   replaced; an exact patch whose face can't be cut to shape is tried as a blend). A
   blend must pass within 0.02 mm of every mesh corner and may not bow away from a facet
   more than a circular arc through its corners would. Patches next to gaps or invalid
   faces left by sewing are dropped individually (blends first, then the nearest patch);
   if the volume check fails, the feature list is halved to find the culprits (for ten
   minutes at most; whatever is left then stays faceted). A patch that runs all the way
   round its axis but is cut to shape is built as two half rings if it won't build whole.
   A final tidy-up is kept only if the solid is still valid. If the STL itself is not a
   clean solid (it touches or crosses itself), the result is checked against what its
   bare facets give and flagged with a warning.

Limits: large freeform areas (organic shells, variable-radius rounds over big areas)
mostly stay faceted: the smooth patch fitter manages small blends and corners, but not
yet whole freeform surfaces. Such parts also convert slowly (tens of minutes at ~17k
triangles). Screw threads (helical surfaces) stay faceted.

## Tests

```
python tests/run_tests.py
```

Builds five CAD test parts (holes, pins, chamfers, crossing holes, a filleted lid,
domes, balls...), meshes them, converts them back and checks each against the original
design (volume, validity, face count).

## Credits

The permissive final pass, the outline-polygon fallback face builder, torus-from-strips
and radius/axis snapping follow ideas
from [stlToSolid](https://github.com/Crypto69/stlToSolid) (no code copied). Built on
[OpenCascade](https://dev.opencascade.org/) via
[cadquery-ocp](https://pypi.org/project/cadquery-ocp/).
