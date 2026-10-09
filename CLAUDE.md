# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

stl2curves converts STL meshes into solid STEP files whose curved areas are rebuilt as
true surfaces (planes, cylinders, cones, spheres, tori, screw threads, smooth freeform
blends), so the part can be edited in Fusion or FreeCAD. It is built on OpenCascade via
`cadquery-ocp`. The owner intends to open-source it (possibly as a FreeCAD plug-in), so
keep dependencies to pip-installable, permissively licensed packages and write methods
from published papers ourselves rather than copying GPL code.

## Commands

Python 3.14, main install, no venv (deliberate: the tool runs permanently from this
clone via a Windows "Send To" shortcut to `stl2curves.bat`, which puts the clone on
`PYTHONPATH` and runs `python -m stl2curves`). The code is the package `stl2curves/`
(`convert.py` is the pipeline and CLI; modules import each other relatively), installable
with `pip install .` / `pip install git+https://github.com/Duvdivan/stl2curves`, which
gives the `stl2curves` and `stl2curves-render` commands. Keep it a package: installed
flat, `build.py`, `features.py`, `workers.py`... would land in site-packages as top-level
modules (`build` clashes with the PyPI `build` tool). The version is
`stl2curves/__init__.py` `__version__` (pyproject reads it); bump it for each release
tag. Supported Python: 3.11-3.14 (cadquery-ocp 8's wheels; it also pulls in vtk).

```
pip install -r requirements.txt                 # cadquery-ocp>=8 (+vtk), numpy, scipy
python -m stl2curves part.stl                   # -> part.step next to it (from the clone root)
python -m stl2curves part.stl --out DIR --details   # list every rebuilt / skipped feature
python -m stl2curves.render part.step           # PNG preview coloured by surface type
python tests/run_tests.py                       # regression suite (5 CAD parts, a few minutes)
```

From another directory, set `PYTHONPATH` to the clone root first (as `stl2curves.bat` and
`compare/regress.sh` do). Scripts that import it use `import stl2curves.convert as S`
(`S.stl_to_solid`, `S.attempt`...); `from stl2curves import stl_to_solid, count, volume,
write_step` also works (loaded lazily, so `import stl2curves` alone doesn't load OCC).

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

- `C:\Users\delta\Downloads\X2D+Accessory+Toolbox_stls` (4 STLs): baseline 114 / 188 /
  211 / 120 faces (s2), all valid. Re-check after every change.
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

`stl2curves/convert.py` (`stl_to_solid`) runs the pipeline (all modules below are in
`stl2curves/`):

1. **`repair.py`** mends the mesh (duplicates, zero-thickness pairs, slivers, cracks,
   orientation, holes, small self-crossings). **`simplify.py`** thins meshes over
   150k triangles with a strict error bound (vertices never move off the surface).
2. **`bodies.py`** splits bodies that only touch along an edge; each is built on its own
   and fused at the end. A 3MF object's parts (`read3mf` labels each triangle) are
   likewise repaired and built one by one and fused (`_repair_parts`): Bambu
   multi-colour objects have parts pressed into pockets of others, and repairing the
   merged mesh cut away their shared walls and left it crossing itself.
3. **`features.py` `analyze()`** groups triangles into flat facets (`Mesh`), then finds
   patches in passes over the smooth regions (`Region`, whose `free` flags track what is
   still unexplained). The pass order matters:
   - Pass 0: whole regions as one surface (strictly; failing that, a cylinder or cone
     with every corner within 0.02 mm, `_whole_noisy`, for coarse exports: cut into
     strips instead, each strip fits a cylinder of its own, wrong radius).
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
5. **`blends.py`** first joins fillets along curved edges (**`pipes.py`**): a
   constant-radius fillet whose edge is neither a line nor a circle (a flat face meeting
   a tilted cylinder: the ball's path is an ellipse) comes out of detection as a chain of
   narrow tori, cylinders, spheres and cones of the fillet's radius, each centred
   elsewhere. A chain grows from its longest piece over same-radius neighbours (each kept
   if one tube still fits all: a tube of radius r round a cubic B-spline spine, fitted by
   point-distance minimisation, each step fitting the spline through the ball centres at
   the corners' foot points; Gauss-Newton on the distances alone let the spine slide
   along itself and blow up) and over loose facets and pieces already on the tube (ones
   reaching past the tube's end join only if it refits with them: the face is cut from
   the tube, which must run on past its outline). A piece is left alone if it is a real fillet:
   the two faces along its sides agree with its axis (`_real`: planes square to a torus's
   axis, cylinders round it; planes and cylinders along a cylinder's). A spline can't
   tell (it bridges a line-to-arc junction within 0.002 mm once the arc is 3 mm or more),
   and reach can't either (round an ellipse's ends a torus fits for 90 deg). Pipes are
   `kind="trimmed"` with a `Pipe` model (`kind="pipe"`, built like freeform faces); one
   that fails to build gives back its pieces.
   Then `blends.py` groups leftover curved facets (and "scrap" pieces that only stand in
   for a surface) into smooth areas. A big area (over `MAX_FACETS` facets or 5 mm
   across) is first tried as one **`freeform.py`** surface: a P-spline height field
   over a plane or a cylinder, fitted to every corner within `FIT_DEV` (0.002 mm), knots
   refined only as needed (at most one control point per corner, spans per direction
   by corner density), bulge-checked between corners like blends. An area that doesn't
   fit is halved by facet direction (`_freeform_pieces`); narrow exact strips between
   freeform pieces are absorbed when one surface fits all (`_absorb_strips`: Fusion's
   curvature-continuous fillets are splines whose middle fits a cylinder over a narrow
   band); loose facets already on a piece's surface are taken in (`_grown`: fringes
   left by halving have edges shorter than any sewing tolerance). Freeform patches are
   `kind="trimmed"` with a `Freeform` model (`signed`, `normal` like any surface); a
   big one blamed for trouble loses only its facets near it (`blends.carve`, trouble
   points kept on `mesh.trouble` by `culprits`), one that fails to build falls back to
   the 60-facet blends of its area (`blends.fallback`). Whatever no freeform takes
   becomes N-sided blends.
6. **`build.py`** makes one face per feature plus a planar face per remaining flat facet,
   sharing each boundary edge between the two faces beside it (`Edges` cache), then
   sews. `convert._build`/`attempt` then checks the solid: valid, and volume equal to
   the faceted volume plus each feature's predicted `change` (within the summed
   `tolerance`; OCC's volume, re-measured by `tessellated_volume` when it disagrees),
   and its STEP file read back the same (`file_trouble`). Features that fail are dropped and left faceted (blends give back the
   pieces they replaced; failed exact patches are retried as blends; whole rings as two
   half rings; a halving search finds culprits, capped in time). After
   `TROUBLE_SECONDS` per body the drop loop "hurries": `culprits` blames every patch
   within 0.5 mm of trouble, four times further each round (a 920k-triangle mesh spent
   95 min dropping a few patches per 3-minute sew). Blame ties (trouble on an edge
   between patches) go to the smallest patch. On a part of 20,000+ faces the final
   tidy-up runs in a worker, given up after `TIDY_SECONDS` (UnifySameDomain on 218k
   faces ran over 20 min); smaller parts are tidied in place, because a shape read back
   from a BRep file in a worker no longer merged at all (Rak N 1780 faces, not 1279).
   A merged face that comes out invalid has its pieces kept apart in a retry
   (`_bad_merges`, `KeepShape` on the edges between them; KeepShape takes edges or
   vertices, not faces), so one bad merge no longer throws away hundreds of good ones.
   A flat facet whose outline touches itself (holes meeting its edge at a corner) is
   built as a few faces over unpinched groups of its triangles (`build._unpinched`), not
   as hundreds of triangles left for the tidy-up to merge.

Speed (`workers.py`): OpenCascade and small-array numpy hold the GIL, so parallel work
runs in a shared process pool (up to 16 workers, fewer if free memory is short: each
holds its own copy of the mesh, ~1.5 KB a triangle, over a gigabyte on a 900k-triangle
mesh), started in the background once a mesh has 3,000+ triangles and kept for further
parts of the same run. Shapes cross processes
as BRep files, the mesh once per pass via `workers.share`/`load`. What runs there:

- The seed passes (pass 1, loose pass, axis passes: `features._seed_pass`), per region,
  and big regions in slices whose answers are replayed in order. The result must stay
  **identical** to the sequential run: check with an analysis fingerprint (md5 of the
  features' sorted facets) on Rak N and the full toolbox after any change there.
- Cutting trimmed faces from their surfaces (`build._cut_all`): outlines are worked out
  in the main process in patch order (`_trimmed_outline`; a run two patches share is
  built once, by whichever comes first), the slow cut (`_trimmed_cut`) in the workers.
  The first build of the GPS case went from 148 s to 12 s.
- Blend fills (each capped at `FILL_SECONDS`) and their fit checks (`_assess_blend`),
  the bare-facet build (`_bare_ahead`), and the final fuse (capped at `FUSE_SECONDS`).
- `TIME_LIMIT` (`--time-limit`) only stops optional refinement (blends, culprit search,
  second chances); analysis always runs to the end (cutting it gave garbage).

FreeCAD add-on (`package.xml` at the root, code in `freecad_addon/`): the repository
itself is the add-on (FreeCAD clones it into `Mod/stl2curves`; `package.xml`'s workbench
`subdirectory` makes FreeCAD run `freecad_addon/InitGui.py`, which registers the
`Stl2Curves_Convert` command and a workbench manipulator putting it after
`Part_ShapeFromMesh` / `Mesh_FromPartShape` and on the Part Tools / Mesh Tools
toolbars). The conversion runs out of process (`QProcess`: `python -m stl2curves`) with
FreeCAD's bundled Python (`bin/python.exe`; `sys.executable` is FreeCAD.exe) and
`PYTHONPATH` = a private `pip install --target` folder (`<UserAppData>/stl2curves/py311`,
cadquery-ocp + numpy + scipy) + the repository root: FreeCAD's own OCC (Part, and
pythonocc in its site-packages) is a different build from cadquery-ocp. So **the package
must keep running on Python 3.11** (FreeCAD 1.0/1.1 bundle 3.11). Meshes are converted
in local coordinates and the solid gets the mesh's global placement. Keep
`package.xml`'s `<version>` equal to `__version__`. A portable FreeCAD 1.1.3 is in
`C:\Users\delta\Tools\FreeCAD`; test the add-on with
`FreeCAD.exe --user-cfg <copy of user.cfg> -M <repo> script.py` (loads it without
installing; a script's prints don't reach the console, so write to a file).

Shared corners and rejoin are on by default (since 2026-10-08). The other experiments
(2026-10-07, after a survey of other tools: see memory prior-art) are each off unless named
in the environment variable `STL2CURVES_TRY` (comma-separated), so each can be judged alone
by full regressions (`compare/regress_try.sh` runs them all):

- `shared`, on by default (`build.Edges`, `SHARED_TOLERANCE`, `SHARED_CROWD`): `build_faces` gives every
  edge the one vertex per mesh point and shares flat-to-flat edges, kept for all of a
  part's attempts, so faces arrive at sewing joined (GPS bare sew 11 -> 4 s, rounds 9 ->
  7 s). Blend fills keep their own edges: one fill on shared edges ground on for 20 min.
  Hazards met: anything that mends shapes in place spreads through shared corners (the
  splitter widened them to 8 mm: `SetNonDestructive`; ShapeFix on flat faces is tried on a
  copy first, `SHARED_TOLERANCE`; the checks after sewing run on a copy). Mesh points within
  `SHARED_CROWD` of another (duplicates, slivers 0.0001 mm across) keep corners of their own
  (`_crowded`): sewing merges such points for unjoined faces, and faces joined there kept
  the tiny edges and no longer met their unshared neighbours (ratchet: no solid). Don't
  close up edges shorter than the sewing tolerance yourself: fine meshes have real edges of
  0.005-0.05 mm (the tensioner came out invalid).
- `rejoin`, on by default (`build._rejoin` in `sew`, `REJOIN_SNAP`): on sew's copy, corners at one mesh
  point (within `REJOIN_SNAP`: copies read back from files) become one vertex, and free
  edges with the same two corners and the same curve become one edge (the curved face's
  kept, the other face's curve moved onto its surface; a flat face takes a line or arc
  either way round). Arc-beside-chord pairs are gaps, left to sewing. Two rules keep it
  from breaking sewing. An edge shorter than twice the sewing tolerance is left to sewing,
  with every edge meeting it: sewing closes the short edge up, and joined round it the edge
  came out closed in one face and gone from the other (no closed solid; the halving search
  took Rak SMA-USB to 1,754 faces). And the corners of a short edge stay apart: merged, they
  changed how sewing matched the edges round them (Uni-USB got a free edge at the end of a
  rounded edge, which the blame then dropped: +50 faces). GPS first round: free edges
  19,650 -> 9,340, sew 8.4 -> 4.5 s.
  Default (d1) vs all off (n7), same day: whole set of 28 parts 14,775 vs 14,753 faces,
  2,359 vs 2,349 s; Rak Uni-USB 1,531 (as all off); Solar Panel Bracket Large 2,256 faces,
  valid (all off: 2,222, invalid); the ratchet is invalid in both (d1: 1,753 faces, with a
  STEP read-back warning that all off doesn't have). GPS 8,242 faces / 771 s vs 8,598 /
  785 s. Before the corner rule (n6): whole set 14,780 faces, 2,286 s; Uni-USB +50 faces;
  bracket 2,215 valid; GPS 8,195 / 675 s. So the corner rule buys back the Uni-USB fix at
  about 73 s over the set (96 s on GPS) and +41 faces on the bracket, netting -5 faces over the set.
- `exact` (build.EXACT_EDGES): a run between two analytic surfaces (a flat facet's plane,
  cylinder, cone, sphere, torus) is their exact intersection (`GeomAPI_IntSS`), trimmed
  at the run's corners (widened to cover their distance off it), if every point of the
  run lies within `EXACT_DEV`. A periodic curve on a single mesh edge must go the short way
  (a 0.2 mm chord became a 1,055 mm ellipse), and the length must match the run's.
- `swap` (assemble.py): no sewing and no rounds. The bare facets' faces on shared edges
  make the start; each patch's face is made on its surface from the edges round its
  facets (ShapeFix adds pcurves; tolerances restored if refused), checked alone (valid,
  area, outward, volume roughly, edges reused, `MAX_EDGE_TOL`) and swapped in, else its
  fallbacks (as `_build` picks them) are tried. Version 0 kept the mesh's chords as the
  curved faces' outlines (GPS: 7,158 faces, valid, read back OK, build stage ~315 s; but
  regression 18,824 faces, turned parts worst). Version 1: whole rings as half rings
  first, the solid checked as `attempt` does and blamed patches taken out again in place
  (`SWAP_ROUNDS`, `Registry.revert`), runs along a flat or swapped face one spline in
  both faces (`SMOOTH_RUNS`), spindle tori swept, surfaces trimmed to the patch. head:
  340 faces vs 40 sewn: a patch is judged before its curved neighbours exist, so single
  chords against a coarse neighbour stay, beyond `MAX_EDGE_TOL`. Every fix here redoes
  what build's face construction already does: the likelier route is build's faces,
  joined (`rejoin`), without the sewing.
- `local.py` (helpers, always available): swap faces/edges in a big shape at the cost of
  the faces touched (`apply_local`: two-step ReShape), local edge maps; `tests/test_local.py`.

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
- `file_trouble` matches faces with a KD-tree (`table`/`partner` in `convert.py`): the
  all-pairs search took 12.7 s on 30,000 faces, the tree 0.45 s, with the same partners.
  The read-back at the end of `_stl_to_solid` can't be skipped: `attempt` checked the shape
  before `_tidy`, and a tidied shape is a new one with a new STEP file.
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
- `repair` cuts out small self-crossing knots and patches the hole; on a fine mesh the
  "knot" can be thousands of triangles whose hole doesn't patch, which left a watertight
  900k-triangle mesh open and unbuildable. A cut is kept only if its patch closes.
- Splines through a run of mesh corners swing far out where the corners are unevenly
  spaced (dense round a bend, one 20 mm step on the straight: 0.57 mm off), and no
  surface then contains the edge; `build._even` adds points along long steps.
- OCC's BRepMesh tessellates a big B-spline face coarsely whatever deflection is asked
  (0.07 mm off a 64x13-span surface at 0.002), so `deviation.py` reports such faces
  worse than they are; measure against the surface itself (GeomAPI_ProjectPointOnSurf).
- `GeomAPI_PointsToBSplineSurface` is unreliable on a sampled grid: allowed degrees up
  to 8 it can report success while missing its own points by millimetres (a sampled
  tube, a cylinder-based freeform surface), and kept cubic it can wave a tenth of a
  millimetre between them (the tray's pockets: 268 -> 360 faces). No one setting is
  right for every surface, so `freeform.approximate` tries a few (chord-length or
  iso-parametric, degree 8 or 3), measures each midway between the grid points, and
  keeps the first within 0.0004 mm.
- A smooth area that fits only with a dense knot grid has a crease or a tight round
  inside it: capping control points by corner count (and splitting) beat refining
  (131x131 control points made the splitter take minutes per face).
- Never hand shared outline edges to ShapeFix (or anything that mends shapes): it widens
  their tolerances in place, even when the face it builds then fails, and the faces
  beside them change with it. That made every build depend on the order faces were
  built in (the same cut twice gave different faces; workers, with fresh copies, gave
  different parts). `_outline_face` works on copies. Check with a cut-twice test.
- OCC's splitter on a degree-8 swept-fillet surface took 78 s to fail where laying the
  outline on it took 0.1 s: pipes get their face from their outline first.
- Sewing costs ~1 s per 1,000 faces whatever the options; faces share no edges before
  it. Incremental sewing (keep the last round's sewn shell, swap in the changed faces,
  local-sew them with `Load`/`Add`, mapping faces through `ModifiedSubShape`, not
  `Modified`) sews a 10,000-face round in ~1 s, but was tried three ways (2026-10-07)
  and lost every time: a sewn face's edges carry the gaps it closed as tolerance, so
  sewn faces join new ones that a fresh sew leaves apart (2 free edges where a full sew
  of the same faces left 41), forced-shut gaps come out invalid, and faces never sewn
  afresh keep their defects round after round. The blame rounds then went differently:
  GPS 8,598 faces / 666 s became 17,099 / 856, 9,163 / 682 (forgetting the shell after
  a failed check) and 8,889 / 702 (re-sewing the faces round the change unsewn). The
  attempt loop is tuned to what a fresh sew does; don't try this again without making
  it give the same free edges as a full sew.
- Blaming every patch near trouble that recurs (85-100% of each round's trouble lies
  within 1 mm of an earlier round's) cut rounds but lost curves and no time (GPS 8,598
  -> 11,805 faces): more faces make every later sew slower.
- `extrude.py` (profiles pushed along a design direction) is switched off
  (`ENABLED`): on the GPS case its few patches mostly failed and took neighbours with
  them. To try on parts made of extruded profiles (gears, cams).
- ShapeFix_Solid can split a zero-volume bubble (two faces lying on each other) off as a
  second shell of the solid: invalid with every face valid and nothing to blame
  (`_without_empty_shells`).
- Big multi-body meshes sew into hundreds of zero-volume shells (coincident face pairs);
  `solids_from_shells` drops them (`EMPTY_SHELL`), or body counts and fuses go wrong.

## Conventions

- Docstrings and comments explain *why* in plain words; module docstrings describe the
  method (with paper references where one is followed). Tunable constants sit at the
  top of each module with a one-line comment giving the unit and meaning.
- Commits: plain-sentence subject lines (e.g. "Find fillets between flat faces by
  rolling-ball fitting"), committed with the repo-local noreply email; remote is the
  public GitHub repo `Duvdivan/stl2curves` (MIT licence), branch `main`: never commit
  downloaded meshes or other people's models.
- Keep `README.md` (user-facing: what is rebuilt, options, how it works, limits) in step
  with behaviour changes.
