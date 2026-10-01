# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

stl2curves converts STL meshes into solid STEP files whose curved areas are rebuilt as
true surfaces (planes, cylinders, cones, spheres, tori, screw threads, smooth freeform
blends), so the part can be edited in Fusion or FreeCAD. It is built on OpenCascade via
`cadquery-ocp`. The owner intends to open-source it (possibly as a FreeCAD plug-in), so
keep dependencies to pip-installable, permissively licensed packages and write methods
from published papers ourselves rather than copying GPL code.

## Commands

Python 3.14, main install, no venv (deliberate: the tool is installed permanently and
runs from a Windows "Send To" shortcut via `stl2curves.bat`).

```
pip install -r requirements.txt                 # cadquery-ocp>=8, numpy, scipy, vtk
python stl2curves.py part.stl                   # -> part.step next to it
python stl2curves.py part.stl --out DIR --details   # list every rebuilt / skipped feature
python render.py part.step                      # PNG preview coloured by surface type
python tests/run_tests.py                       # regression suite (5 CAD parts, a few minutes)
```

`tests/run_tests.py` generates its STLs into `tests/parts/` (gitignored) with the
`make_*.py` scripts, converts them and checks validity, volume error and a maximum face
count. There is no per-test selector; to run one case after the parts exist:

```
python -c "import sys; sys.path.insert(0,'.'); from stl2curves import stl_to_solid, count, volume; from OCP.TopAbs import TopAbs_FACE; s,_=stl_to_solid('tests/parts/test_lid.stl',0.01,fuse=False); print(volume(s), count(s,TopAbs_FACE))"
```

There is no linter or formatter configured.

## Judging a change

The unit suite is not enough. Detection changes routinely look harmless on radii but
cost 30-50% more faces downstream, so judge every change by **full conversions** of the
real-world parts and compare face count, validity and time against the last run. Never
edit code while a batch of conversions is running: processes started later import the
edited code.

Local regression data on this machine (not in the repo):

- `C:\Users\delta\Downloads\X2D+Accessory+Toolbox_stls` (4 STLs): baseline 203 / 194 /
  243 / 230 faces, all valid. Re-check after every change.
- `C:\Users\delta\Tools\compare\rak` (20 parts, coordinates rounded to 0.001 mm):
  `run2.sh`, `summary.sh`, `devall.sh` there; outputs in `out_v*`.
- `compare\bass` (instrument: bent tubes, cavities), `compare\box2` (`tray.stl`, and
  `toolbox.stl`: 223k triangles, many bodies touching at corners, ~11 min; both with the
  original STEPs as ground truth), `compare\new5` (newer downloads, incl. a Bambu 3MF).
- `C:\Users\delta\Tools\compare\regress.sh`: Rak, toolbox parts, instrument and tray,
  one conversion at a time (each uses all the worker processes; running several at
  once only makes the timings meaningless).
- `C:\Users\delta\Tools\compare\gps` (freeform-heavy, 20+ min per part) and
  `compare\repeater` (1.3M-triangle hubs, exercises `simplify.py`): slow; only run
  them when the change targets them.
- Mesh-to-solid deviation: `C:\Users\delta\Tools\compare\stlToSolid\.venv\Scripts\python.exe C:\Users\delta\Tools\compare\deviation.py part.stl part.step`

## Architecture

`stl2curves.py` (`stl_to_solid`) runs the pipeline:

1. **`repair.py`** mends the mesh (duplicates, zero-thickness pairs, slivers, cracks,
   orientation, holes, small self-crossings). **`simplify.py`** thins meshes over
   150k triangles with a strict error bound (vertices never move off the surface).
2. **`bodies.py`** splits bodies that only touch along an edge; each is built on its own
   and fused at the end.
3. **`features.py` `analyze()`** groups triangles into flat facets (`Mesh`), then finds
   patches in passes over the smooth regions (`Region`, whose `free` flags track what is
   still unexplained). The pass order matters:
   - Pass 0: whole regions as one surface.
   - Screw threads (`threads.py`), before anything tries surfaces round known axes or
     the loose pass (which would take thread pieces for cylinders and cones). A thread
     must reach at least 270 deg round its axis and be no deeper than its lead: a finely
     meshed smooth bend agrees with some screw motion over a narrow arc (a 58k-triangle
     bracket gave three false threads and spent minutes on 130 candidates).
   - Pass 0b: surfaces on the axes found so far (a lug round its screw hole).
   - Fillets between flat faces (`fillets.py`, rolling ball): before pass 1, which
     would cut a coarsely meshed fillet into strips of its own. Forgiving (most band
     corners within 0.02 mm, all within 0.05), so its guards matter: faces at least 3x
     the band's facets (finely cut pins), faces and fillet running at least one radius
     along the edge (slices of rounded corners and curved edges), and radii found with
     confidence preferred for near-misses (one radius is used all over a design).
     Everything is done in a cross-section (`_Flat`, `_Wall`, `_Beside`) with the ball
     rolling between two curves (`_Lines`, `_LineCircle`). After flat-flat pairs come
     flat face + wall pairs, the walls being the features so far (pieces on one
     surface joined): face square to the axis gives a torus, face along a cylinder's
     axis a cylinder. A coarse chamfer round a hole fits a torus exactly (all its
     corners on the rims), hence the profile-bend and interior-corner guards.
   - Pass 1: surfaces hypothesised from pairs of neighbouring facets.
   - Pass 2: tori around axes already found; thread ends (countersinks, chamfers).
   - Loose pass, band tori, merge of same-surface neighbours.

   Thread features are appended after the merge: their models aren't surfaces of
   revolution.
4. **`sizing.py`** guesses the design size (one scale factor making radii round);
   `features.snap()` equalises radii and snaps near-axis-aligned axes.
5. **`blends.py`** groups leftover curved facets (and "scrap" pieces that only stand in
   for a surface) into smooth freeform blend regions.
6. **`build.py`** makes one face per feature plus a planar face per remaining flat facet,
   sharing each boundary edge between the two faces beside it (`Edges` cache), then
   sews. `stl2curves._build`/`attempt` then checks the solid: valid, and volume equal to
   the faceted volume plus each feature's predicted `change` (within the summed
   `tolerance`; OCC's volume, re-measured by `tessellated_volume` when it disagrees),
   and its STEP file read back the same (`file_trouble`). Features that fail are dropped and left faceted (blends give back the
   pieces they replaced; failed exact patches are retried as blends; whole rings as two
   half rings; a halving search finds culprits, capped in time).

Speed (`workers.py`): OpenCascade and small-array numpy hold the GIL, so parallel work
runs in a shared process pool (up to 16 workers), started in the background once a mesh
has 3,000+ triangles and kept for further parts of the same run. Shapes cross processes
as BRep files, the mesh once per pass via `workers.share`/`load`. What runs there:

- The seed passes (pass 1, loose pass, axis passes: `features._seed_pass`), per region,
  and big regions in slices whose answers are replayed in order. The result must stay
  **identical** to the sequential run: check with an analysis fingerprint (md5 of the
  features' sorted facets) on Rak N and the full toolbox after any change there.
- Blend fills (each capped at `FILL_SECONDS`) and their fit checks (`_assess_blend`),
  the bare-facet build (`_bare_ahead`), and the final fuse (capped at `FUSE_SECONDS`).
- `TIME_LIMIT` (`--time-limit`) only stops optional refinement (blends, culprit search,
  second chances); analysis always runs to the end (cutting it gave garbage).

Key contracts:

- **`Feature`** (dataclass in `features.py`): `model`, `kind` (`revolve`, `wedge`,
  `ball`, `trimmed`, `blend`), `facets`, and `change`/`tolerance`/`worst`, which feed the
  volume check and the sewing tolerance.
- **Surface models** (`Revolved`, `Sphere`, `threads.Helical`/`threads.Flank`) expose
  `signed(p)` (distance, sign fixed per model), `normal(p)` consistent with that sign,
  and `kind`. Code that assumes an axis or radius must guard with
  `isinstance(m, Revolved)` / `Sphere`.
- **`trimmed` faces** are built by splitting a generous piece of the surface
  (`build._generous_surface`) with the patch's outline edges and keeping the pieces the
  facets' projections land on; `Boundary.labels` marks outline runs that follow a
  surface's natural boundary lines, so they get exact lines and arcs instead of splines.
- Tolerances are per mesh (`features._mesh_tol` = 1 µm + 2e-7 x largest coordinate).
  Files that round coordinates (often to 0.001 mm) are detected by `grid_noise`, and
  that noise is allowed for when merging coplanar triangles and in surface tests.

## Hard-won lessons (don't undo these)

- The boolean "skin" approach was fragile; faces are built and sewn instead. Don't go
  back to it.
- Two-facet cylinder axis guesses fail on irregularly triangulated fillets; the
  least-squares fallbacks and permissive-pass scraps exist for a reason (removing them
  costs faces).
- Blend acceptance: within 0.02 mm of mesh corners, bow no more than a circular arc
  through its corners (`blends.py` constants). Guide points made blend surfaces wave, so
  they are off. A blend's sewing gap is its edge tolerance, not its mesh deviation.
- Widening detection tolerance by the grid noise produced false spheres; so did global
  noise-tolerant sliver handling (broke the RAK enclosures).
- Facet normals on real meshes are good only to about 1°. Anything needing precision
  (thread pitch, axis) should come from corners, which are exact to the file's rounding.
- OCC's volume integration (GProp) goes astray on solids sewn across wide gaps (a
  curved patch's outline edge bowed onto its surface, the flat facet beside it keeping
  the chord, 0.05-0.09 mm apart on a rough mesh): several mm^3 off, the same wherever
  the part sits, while a watertight tessellation gives the predicted change. That is
  why `attempt` re-measures by tessellation before failing a volume check.
- An in-memory solid that checks out can still write a broken STEP file (12 of 28
  regression parts did, 2026-10-01): sewing and ShapeFix widen edge tolerances up to
  ~0.3 mm, which hides edges crossing inside a face; STEP stores no tolerances and the
  reader recomputes tight ones; and an edge shorter than the sewing tolerance gets both
  ends merged, which a reader takes for a whole circle. Hence `file_trouble` in
  `attempt` (write, read back, compare) and `build._without_collapsed` after sewing.
  Check outputs with a read-back, not only `BRepCheck_Analyzer` on the shape in memory.
- OCP 8 quirks: `TopoDS.Shell/Face` (no `_s`), `OCP.collections` for arrays and
  sequences, `Bnd_Box.Get()` is broken, `Quantity_Color` returns linear RGB,
  `BRepCheck_Result.Status()` can't be read (use `BRepCheck_Analyzer.IsValid(sub)`).
- Pickling a `Mesh` reorders its neighbour *sets*, and pass 1 tries `nbrs[i][:3]`: send
  regions' neighbour lists along to workers, or results silently change.
- Windows process pools spawn a worker only when a job finds none idle, one by one in
  the submitting thread: `workers._start` calls `_launch_processes()` up front. Scripts
  that convert must have an `if __name__ == "__main__"` guard (workers re-import them).
- Profiling: cProfile inflates small functions; py-spy on Windows hangs or samples the
  launcher. A sampling thread over `sys._current_frames()`, weighted by elapsed time
  (long OCC calls hold the GIL), works well.
- Big multi-body meshes sew into hundreds of zero-volume shells (coincident face pairs);
  `solids_from_shells` drops them (`EMPTY_SHELL`), or body counts and fuses go wrong.

## Conventions

- Docstrings and comments explain *why* in plain words; module docstrings describe the
  method (with paper references where one is followed). Tunable constants sit at the
  top of each module with a one-line comment giving the unit and meaning.
- Commits: plain-sentence subject lines (e.g. "Find fillets between flat faces by
  rolling-ball fitting"), committed with the repo-local noreply email; remote is the
  private GitHub repo `Duvdivan/stl2curves`, branch `main`.
- Keep `README.md` (user-facing: what is rebuilt, options, how it works, limits) in step
  with behaviour changes.
