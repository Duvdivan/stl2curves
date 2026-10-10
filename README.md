# stl2curves

Convert STL meshes and 3MF projects (3D-printing downloads, slicer projects, exports
from other programs) into solid STEP files with **true curved surfaces**, so they can
be edited properly in Fusion, FreeCAD or any other CAD program.

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
| swept fillets | constant-radius rounded edges along any other curve (where a flat face meets a tilted hole or slot, along a spline-shaped outline) |
| screw threads | bolts, threaded holes and sockets (helical surfaces, with the thread named) |
| freeform surfaces | crowns that fade out round a bend, curvature-continuous ("smooth") fillets, walls of slots milled along a curve |
| smooth blends | where rounded edges meet at a corner, fillets along spline-shaped edges |

Surfaces can be cut to any outline (a hole running out through a sloped face, two
holes crossing), and overlapping bodies in one STL are fused automatically, as are
bodies that only touch along an edge.

Curved areas that none of the exact surfaces fit become freeform faces: a big smooth
area one B-spline surface fitted to all its mesh corners (within 0.002 mm), a small one
(where rounded edges meet at a corner) one smooth patch spanning its outline. Anything
else stays as flat facets, so the result always matches the mesh.

> **Beta.** stl2curves works on the parts it has been tested on (mostly 3D-printing
> downloads and CAD exports), but every mesh is different. If a conversion fails, takes
> very long or loses curves it should have kept, please
> [report it](https://github.com/Duvdivan/stl2curves/issues/new/choose).

## Install

You need Python 3.11, 3.12, 3.13 or 3.14 (from [python.org](https://www.python.org/downloads/);
on Windows tick "Add python.exe to PATH" in the installer). Then, in a terminal:

```
pip install git+https://github.com/Duvdivan/stl2curves
```

(Without git: download the repository as a ZIP from GitHub, unzip it and run
`pip install .` in that folder.) This installs the `stl2curves` command and its
dependencies: OpenCascade (`cadquery-ocp`, which brings VTK along), numpy and scipy,
about 750 MB in all.

Update with `pip install --upgrade --force-reinstall git+https://github.com/Duvdivan/stl2curves`,
and check the version with `stl2curves --version`.

## Use

```
stl2curves part.stl                    # -> part.step next to it
stl2curves project.3mf                 # one STEP per object printed in it
stl2curves some_folder --out STEP      # every .stl and .3mf in a folder
stl2curves a.stl b.stl --merge all     # also one combined STEP
stl2curves part.stl --details          # list every rebuilt feature
stl2curves part.stl --true-size        # rebuild at the apparent design size
stl2curves-render part.step            # PNG preview coloured by surface type
```

`python -m stl2curves ...` does the same as `stl2curves ...` (handy if the command isn't
on your PATH).

On Windows, `stl2curves.bat` (in the repository) accepts drag-and-drop of STL/3MF files
or folders onto it; put a shortcut to it in `shell:sendto` for a right-click "Send to"
entry. It runs the copy of stl2curves next to it, so it also works from an unzipped
download (after `pip install -r requirements.txt` there) without installing the package.

Options: `--tol` sewing tolerance (mm), `--no-fuse` keep overlapping bodies separate,
`--no-curves` flat faces only, `--no-blends` no smooth freeform faces (exact surfaces
and flat facets only), `--no-repair` take the mesh as it is, `--simplify MM` thin out
the mesh first (moving its surface by at most MM; meshes over 150,000 triangles are
thinned by 0.005 mm automatically), `--time-limit SECONDS` (default 600, 0 for none),
`--accuracy MM` how far smooth surfaces, fillets on rough meshes and other fitted surfaces
may stray from the mesh (default 0.02): a looser setting rebuilds more of the part as
curves instead of leaving it as triangles. On the GPS case back used in testing, 0.1
gave 857 faces instead of 1,502, with the mesh still within 0.016 mm of the solid at
99% of its corners (0.044 mm at most).

From a 3MF file (Bambu Studio, OrcaSlicer, PrusaSlicer...) every object on the build
plates is converted to its own STEP, named after the file and the object, placed as on
the plate. Modifier volumes, negative volumes and support blockers are left out. An
object made of several parts (a multi-colour print, say) has each part mended and
rebuilt on its own, and the parts joined into one solid at the end (they may touch or
overlap: merged into one mesh first, they would make a mesh that crosses itself).

### Studio: check and rescue areas by hand

```
stl2curves-studio part.stl             # or: python -m stl2curves.studio part.stl
```

opens a page in your browser (served from your own computer only) showing the mesh
coloured by what was found on it: holes, rounded edges, chamfers, rounded corners,
smooth blends. **Convert** builds the solid and shows its faces; areas that came out as
triangles (a fillet it missed, say) are then shown red on the mesh. Paint such an area
with the brush (or pick whole features), one group per face it should become, and
Convert again: each group is rebuilt as one cylinder, cone or sphere if it lies on one,
else as one smooth surface. **Save STEP** downloads the result. The first stages of the
conversion are kept between runs, so converting again after painting skips the
analysis. The page draws with three.js, fetched from a CDN the first time.

### In FreeCAD

stl2curves is also a FreeCAD add-on (FreeCAD 1.0 or later). Until it is listed in the
Addon Manager, install it by hand: download the repository (green "Code" button,
"Download ZIP"), unzip it into FreeCAD's `Mod` folder (in FreeCAD: Macro > Macros...
shows the user macro folder; `Mod` sits next to it, e.g.
`%APPDATA%\FreeCAD\v1-1\Mod` on Windows) so that you get `Mod/stl2curves/package.xml`
(rename the unzipped `stl2curves-main` folder to `stl2curves`), and restart FreeCAD.

Select one or more mesh objects and choose **Part > Mesh to Curved Solid (stl2curves)**
(or **Meshes > ...** in the Mesh workbench, or its toolbar button); with no mesh
selected it asks for STL or 3MF files instead. The solid is added next to the mesh, at
its placement, and the mesh hidden. The conversion runs in the background with its
output shown, and Stop ends it.

On first use it offers to download what it needs (OpenCascade as `cadquery-ocp`, numpy
and scipy, about 750 MB on disk) into a folder of its own in FreeCAD's user data
(`stl2curves/py311` for FreeCAD's Python 3.11); FreeCAD's own packages are not changed.
If you already have a Python with stl2curves' dependencies (3.11 or later), you can
choose it under Options instead.

### Time

Small parts take seconds, typical printed parts of 5,000-15,000 triangles a minute or
two, and very large ones (hundreds of thousands of triangles) about ten minutes. The
slow steps run in parallel on all but two of the computer's cores (up to 16 worker
processes). After the time limit, fitting further freeform faces and hunting down
patches that spoil the solid stop, and whatever has checked out by then is kept: the
result is always a solid matching the mesh, with fewer curves rebuilt the earlier it
had to stop.

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

1. **Facets** — the mesh is mended first (`repair.py`: duplicate and zero-thickness
   triangles, slivers, cracks, flipped triangles, holes and small self-crossings), and
   meshes over 150,000 triangles are thinned out (`simplify.py`) without moving their
   surface more than 0.005 mm. Triangles are then grouped into flat facets (coplanar
   pieces): sliver triangles (three corners in a line) are split away, and when the
   file rounded its coordinates (many exports write 0.001 mm steps) the rounding noise is
   allowed for, both when deciding which triangles are coplanar (long thin triangles
   tilt noticeably) and in every surface test. Pieces of the mesh that merely touch
   (stacked blocks) are converted separately and fused.
2. **Detection** (`features.py`) — first, each smooth area bounded by sharp edges is
   tried as a whole (a countersink, a plain hole): one cylinder, cone or sphere fitted
   through all its corners at once, allowing a few outline corners that a neighbouring
   face's triangulation put slightly off the curve. A mesh exported coarsely keeps its
   cylinders only to a hundredth of a millimetre or so; a whole area that turns far
   enough and has every corner within 0.02 mm of one cylinder or cone (nearly all within
   0.01) is taken as that surface rather than cut into strips of wrong radii. Then neighbouring facets propose a
   surface (cylinder, cone, sphere; tori on axes already found), which is grown while
   mesh corners stay within a micron of it. Fillets between flat faces are looked for
   early (`fillets.py`: two faces meeting at an angle, rounded off, give a cylinder
   touching both, so only the radius is unknown and every corner of the rounded band
   pins it down). This finds fillets however coarsely or irregularly they were cut into
   triangles, and forgives meshes that lost some accuracy there (most corners within
   0.02 mm); and since a design usually uses one fillet radius throughout, a band that
   nearly fits a radius found cleanly elsewhere in the part takes that radius. The same
   rolling ball then looks for fillets between a flat face and a curved wall found so
   far: a floor meeting a boss or hole (a torus round its axis), or a wall with a boss
   or groove running along it (a cylinder beside the wall's). Strict
   passes follow (a patch must follow its surface's natural boundary lines, or meet its
   neighbours at creases), then a permissive pass for what is left, then neighbouring
   patches on the same surface are
   merged, chains of short
   cylinder strips or sphere bits are replaced by the torus they approximate, and equal
   radii and near-axis-aligned axes are snapped to exact values (each snap kept only if
   the patch still lies on the mesh).
3. **Swept fillets** (`pipes.py`) — a constant-radius fillet along an edge that is
   neither straight nor circular is the tube a ball sweeps rolling along it, and
   detection covers it with narrow pieces of tori, cylinders and spheres of its radius,
   each fitting a short stretch. Such a chain becomes one tube round a fitted spline
   (the path of the ball's centre), passing within 0.002 mm of every mesh corner. Pieces
   that are real fillets (lying between faces that share their axis, like a torus round
   a hole's rim) stay as they are.
4. **Freeform surfaces and blends** (`freeform.py`, `blends.py`) — curved areas left
   over (unrecognised facets, and small pieces that only stand in for a surface they
   don't really fit, such as narrow cylinder strips across a spline-shaped bend) are
   grouped into smooth areas, not reaching across creases. A big one gets one B-spline
   surface: a height field over a plane (or, for an area curling further round, a
   cylinder) fitted to its mesh corners by penalised least squares (P-splines, Eilers
   and Marx 1996), with knots refined only until it passes within 0.002 mm of every
   corner and doesn't bow between them; an area no surface fits is halved by facet
   direction and each half tried again, and a narrow exact strip lying between freeform
   pieces (the middle of a curvature-continuous fillet fits a cylinder over a narrow
   band) joins them if one surface fits all; loose facets next to a freeform surface
   that already lie on it are taken in. Smaller areas, and what no surface fits,
   are grouped into blend regions that each turn no more than one smooth patch can
   follow.
5. **Build** (`build.py`) — one exact face per patch, one freeform surface cut to its
   area's outline, one smooth N-sided patch per
   blend (through the mesh corners; blends are built on their own first, a blend that
   won't fit is cut in two and retried) and one planar face per flat area, with shared
   edges (exact lines/arcs where a patch's boundary
   follows its natural lines, splines through the mesh points otherwise, split at sharp
   corners, with points added along long straight stretches so the spline doesn't swing
   out between them), sewn into a solid. Bodies that touch along an edge (`bodies.py`) are built
   separately and fused.
6. **Checks** — every result must be a valid closed solid whose volume matches the
   mesh plus the predicted change from the curves (measured again from a fine
   tessellation where OpenCascade's own volume disagrees: its integration goes astray
   where the sewing had to close gaps of a few hundredths of a millimetre, as beside a
   curved patch on a rough mesh); features that fail are left
   faceted instead of spoiling the part (a blend that fails gives back the pieces it
   replaced; an exact patch whose face can't be cut to shape is tried as a blend). A
   blend must pass within 0.02 mm of every mesh corner and may not bow away from a facet
   more than a circular arc through its corners would. Patches next to gaps or invalid
   faces left by sewing, or next to a face that came out flipped, are dropped
   individually (blends first, then the nearest patch), and once the part builds they
   get a second chance together; if the volume check fails, the feature list is halved
   to find the culprits (until the time limit; whatever is left then stays faceted).
   Each try builds and sews the whole part, minutes on a mesh of a quarter million
   triangles, so after ten minutes of this every patch near a problem is dropped at
   once, reaching four times further each round (0.5 mm, 2 mm, 8 mm...): a few more
   tries at most, losing curves only round the trouble, rather than hours of them. A
   big freeform surface blamed for a problem at its edge gives up only its facets near
   the problem; one that fails outright gives way to the smaller blends its area would
   have had. A
   patch that runs all the way
   round its axis but is cut to shape is built as two half rings if it won't build whole.
   A final tidy-up merges faces on one surface (any merge that comes out broken is
   undone and the rest kept) and is kept only if the solid is still valid. If the STL itself is not a
   clean solid (it touches or crosses itself), the result is checked against what its
   bare facets give and flagged with a warning.
   Finally the solid must survive its own STEP file: it is written out and read back,
   and must come back valid with the same volume (except on parts of over 50,000
   faces, where that would take longer than the whole conversion). A STEP file keeps no
   tolerances, so a reader works them out again from the geometry, tightly; a face that
   only checked out within the generous tolerance sewing gave its edges (two edges
   crossing inside it, or a sliver's edge shorter than that tolerance closed up into a
   loop) would come back broken. The patches at the outline that changed are dropped
   like any other troublemaker, for a few rounds at most; whatever they don't cure is
   reported as a warning. (Edges that sewing closes up are also removed straight after
   sewing.)

Limits: a freeform surface is a height field over a plane or a cylinder, so an area
that curls round in two directions at once (an organic shell, a knob) is cut into
several, and their seams follow mesh edges. Lettering and other shapes extruded from
free curves stay faceted.

## Reporting problems

Open an [issue](https://github.com/Duvdivan/stl2curves/issues/new/choose) with the
command you ran, `stl2curves --version`, your operating system, and the text the
conversion printed (run it with `--details` if you can). A mesh to reproduce the problem
helps most, but only attach one you have the right to share (your own design, or one
whose licence allows it): the issue tracker is public.

## Tests

From a clone of the repository:

```
pip install -e .
python tests/run_tests.py
```

Builds five CAD test parts (holes, pins, chamfers, crossing holes, a filleted lid,
domes, balls...), meshes them, converts them back and checks each against the original
design (volume, validity, face count).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for
full text.

The MIT license allows commercial use, modification, distribution, and private use,
while requiring that the copyright notice and permission notice remain in all copies or
substantial portions of the software.

## Credits

The permissive final pass, the outline-polygon fallback face builder, torus-from-strips
and radius/axis snapping follow ideas
from [stlToSolid](https://github.com/Crypto69/stlToSolid) (no code copied). Built on
[OpenCascade](https://dev.opencascade.org/) via
[cadquery-ocp](https://pypi.org/project/cadquery-ocp/).
