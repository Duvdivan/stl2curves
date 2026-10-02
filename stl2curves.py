"""
stl2curves - convert STL mesh files into solid STEP files for Fusion / any CAD.

What it does:
  1. Reads the STL triangles and groups them into flat facets.
  2. Finds curved areas and works out the exact surface each lies on: holes, pins,
     rounded edges (fillets), chamfers/countersinks around round edges, cones,
     domes, dimples, balls and rounded corners.
  3. Builds a real solid from exact faces: one flat face per flat region (a box
     wall is one face you can select, offset, sketch on, press/pull...) and one
     true curved face per curved area, so Fusion sees real circular edges you can
     select, measure, dimension and fillet.
  4. Writes a .step file next to the STL (or into --out folder).

Straight chamfers and flat faces are exact already. Curved areas that don't match
one of the recognised shapes cleanly (freeform surfaces, variable-radius blends,
odd corner blends) stay as small flat facets.

Usage:
  python stl2curves.py file1.stl [file2.stl ...] [--out FOLDER] [--merge NAME]
  python stl2curves.py some_folder            (converts every .stl in it)

  --merge NAME     also write all inputs together into one NAME.step
  --tol X          sewing tolerance in mm (default 0.01)
  --no-fuse        keep overlapping bodies within one STL as separate bodies
  --no-curves      skip curve detection (leave everything faceted)
  --true-size      rebuild at the size the part appears to have been designed at
                   (undoing e.g. a 99% slicer scale or an inch/cm export), with radii
                   snapped to round values
  --details        list every rebuilt feature with its size
  --simplify MM    thin out an over-dense mesh first, moving its surface by at most MM
                   (0: never; by default meshes over 150,000 triangles, at 0.005 mm)
  --no-repair      don't mend mesh defects first
"""
import argparse
import math
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from OCP.TopoDS import TopoDS, TopoDS_Compound, TopoDS_Iterator, TopoDS_Shape
from OCP.BRepTools import BRepTools, BRepTools_WireExplorer
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE, TopAbs_IN,
                        TopAbs_FORWARD, TopAbs_REVERSED)
from OCP.STEPControl import STEPControl_Writer, STEPControl_Reader, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.IFSelect import IFSelect_RetDone
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopLoc import TopLoc_Location
from OCP.collections import List_TopoDS_Shape, IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher

import features as features_mod
from features import load_stl, analyze, summarize, snap, half_rings, Mesh, TOL, time_left
from sizing import guess_size
from build import build_faces, sew, features_near, point_facet_distance, settle_blends
import workers
from bodies import split_bodies
from blends import add_blends, split as split_blend, _blend as as_blend
from repair import repair
from simplify import simplify
from read3mf import read_3mf

SEARCH_SECONDS = 600    # time allowed for hunting down patches that spoil the solid
TIME_LIMIT = 600        # seconds after which the optional refinements stop (--time-limit); see stl_to_solid
FIX_PRECISIONS = (1e-5, 1e-4)   # mm: precisions ShapeFix tries on an invalid solid after its default
EMPTY_SHELL = 1e-6      # closed shells with less volume than this share of the biggest are dropped
FUSE_SECONDS = 30       # bodies whose fuse takes longer than this are handed over side by side
BARE_AHEAD = 5000       # facets from which the bare-facet build is made ahead, in a worker process
WORKERS_FROM = 5000     # triangles from which a part is worth starting worker processes for
BARE_DRIFT = 1e-4       # share of a sound mesh's volume its bare facets may differ by and still be a reference
TESS_DEFLECTION = 0.0005    # mm: how closely a solid's faces are tessellated to measure its volume
FILE_VOLUME = 1e-3      # share of the volume a finished part's STEP file may read back off by
FILE_MATCH = 0.05       # mm (and share): how far a face's middle or outline may move read back from the file
FILE_AREA = 0.5         # mm^2 (plus half a percent): how much a face's area may change read back from the file
FILE_ROUNDS = 3         # rounds of patches dropped for what their STEP file does, at most, per body
                            # where OCC's integration and the expected volume disagree
UNIFY_TOL = 0.001       # mm: faces this close to one surface are merged at the end ...
UNIFY_DEG = 0.1         # ... if their normals agree this closely
UNIFY_VOLUME = 1e-3     # and kept if the volume stays this close (the integration itself
                        # wobbles by a few mm^3 on a big part)
AUTO_SIMPLIFY = 150000  # meshes with more triangles than this are thinned out first ...
SIMPLIFY_ERROR = 0.005  # ... moving the surface by at most this (mm)


def count(shape, kind):
    n, ex = 0, TopExp_Explorer(shape, kind)
    while ex.More():
        n += 1
        ex.Next()
    return n


def volume(shape):
    # adaptive integration: the default's fixed sample points undercount long spline
    # faces (a thread flank winding five turns came out ~80 mm^3 short)
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-6)
    return props.Mass()


def tessellated_volume(shape):
    """The volume a fine tessellation of the shape's faces encloses, or None if a face
    won't tessellate. OCC integrates face by face, each over its own surface; where the
    sewing joined edges a few hundredths of a millimetre apart (a curved patch's outline
    bowed onto its surface, the flat facet beside it keeping the straight chord) the
    faces don't quite meet, and the integration can come out several mm^3 off, one way
    or the other. The tessellation follows the shared edges, so it closes up."""
    BRepTools.Clean_s(shape)
    BRepMesh_IncrementalMesh(shape, TESS_DEFLECTION, False, 0.1, True)
    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    p, q = box.CornerMin(), box.CornerMax()
    centre = np.array([p.X() + q.X(), p.Y() + q.Y(), p.Z() + q.Z()]) / 2     # (less to cancel)
    total = 0.0
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        face = TopoDS.Face(ex.Current())
        ex.Next()
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is None:
            BRepTools.Clean_s(shape)
            return None
        move = loc.Transformation()
        P = np.array([(n.X(), n.Y(), n.Z()) for n in (tri.Node(i).Transformed(move) for i in range(1, tri.NbNodes() + 1))])
        T = np.array([tri.Triangle(i).Get() for i in range(1, tri.NbTriangles() + 1)]) - 1
        P = P - centre
        v = np.einsum("ij,ij->i", P[T[:, 0]], np.cross(P[T[:, 1]], P[T[:, 2]])).sum() / 6
        total += -v if face.Orientation() == TopAbs_REVERSED else v
    BRepTools.Clean_s(shape)
    return float(total)


def fixed(solid):
    fix = ShapeFix_Solid(solid)
    fix.Perform()
    return fix.Solid()


def compound(shapes):
    comp = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(comp)
    for s in shapes:
        b.Add(comp, s)
    return comp


def boolean(op, args, tools):
    lst = lambda xs: (l := List_TopoDS_Shape(), [l.Append(x) for x in xs])[0]
    algo = op()
    algo.SetArguments(lst(args))
    algo.SetTools(lst(tools))
    algo.SetFuzzyValue(1e-5)
    algo.SetRunParallel(True)
    algo.SetUseOBB(True)
    algo.Build()
    if not algo.IsDone():
        raise RuntimeError("boolean operation failed")
    return algo.Shape()


def mesh_volume(mesh):
    a, b, c = (mesh.pts[mesh.tris[:, k]] for k in range(3))
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6)


def solids_from_shells(sewn, inward=True):
    """Turn sewn shells into solids. A shell inside another is a cavity, the rest are bodies
    (inward False: the mesh has no inward-facing shell, so every shell is a body).

    Returns (shape, number of bodies, number of cavities, total volume of material).
    """
    shells = []
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    while ex.More():
        shell = TopoDS.Shell(ex.Current())
        solid = fixed(BRepBuilderAPI_MakeSolid(shell).Solid())   # oriented outwards
        if solid.ShapeType() != TopAbs_SOLID:
            ex2 = TopExp_Explorer(solid, TopAbs_SOLID)
            if not ex2.More():
                ex.Next()
                continue                    # not a closed piece: ignore it
            solid = TopoDS.Solid(ex2.Current())
        v = volume(solid)
        if v < 0:
            solid = TopoDS.Solid(solid.Reversed())
            v = -v
        shells.append((shell, solid, v))
        ex.Next()
    if not shells:
        raise RuntimeError("no closed body found")
    # (a closed shell with next to no volume is two faces lying on each other, where the
    # mesh had a wall of no thickness: no material, and not a body)
    biggest = max(abs(v) for _, _, v in shells)
    shells = [x for x in shells if abs(x[2]) > EMPTY_SHELL * biggest]
    shells.sort(key=lambda x: -x[2])        # biggest first: an enclosing shell comes first

    def points_of(shape, most=24):
        # corners and face centres (a rib spanning wall to wall has every corner buried
        # in the walls, but the middle of its top face is out in the open)
        faces, ex = [], TopExp_Explorer(shape, TopAbs_FACE)
        while ex.More():
            faces.append(ex.Current())
            ex.Next()
        corners, ex = [], TopExp_Explorer(shape, TopAbs_VERTEX)
        while ex.More():
            corners.append(ex.Current())
            ex.Next()
        pts = []
        for k in range(0, len(faces) + len(corners), max(1, (len(faces) + len(corners)) // most)):
            if k < len(faces):
                # (only the faces sampled: working out every face's centre is slow)
                props = GProp_GProps()
                BRepGProp.SurfaceProperties_s(faces[k], props)
                pts.append(props.CentreOfMass())
            else:
                pts.append(BRep_Tool.Pnt_s(TopoDS.Vertex(corners[k - len(faces)])))
        return pts

    def inside(body, pts):
        # (one classifier per body, set up once: making one per point re-reads a big
        # body's faces each time; a point outside its box is outside it)
        lo, hi = body[2]
        if any(not (lo[0] <= p.X() <= hi[0] and lo[1] <= p.Y() <= hi[1] and lo[2] <= p.Z() <= hi[2])
               for p in pts):
            return False
        if body[3] is None:
            body[3] = BRepClass3d_SolidClassifier(body[0])
        for p in pts:
            body[3].Perform(p, 1e-6)
            if body[3].State() != TopAbs_IN:
                return False
        return True

    def box(shape):
        b = Bnd_Box()
        BRepBndLib.Add_s(shape, b)
        p, q = b.CornerMin(), b.CornerMax()
        return (p.X(), p.Y(), p.Z()), (q.X(), q.Y(), q.Z())

    # A cavity lies wholly inside its host; an overlapping body (a rib sunk into a
    # floor, say) has some corners outside, and gets fused on instead.
    bodies, cavities, total = [], [], 0.0      # bodies: [outer solid, [cavity solids], box, classifier]
    for shell, solid, v in shells:
        pts = points_of(shell) if bodies and inward else None
        host = next((b for b in bodies if inside(b, pts)), None) if pts else None
        if host is None:
            bodies.append([solid, [], box(solid), None])
            total += v
        else:
            host[1].append(solid)
            cavities.append(solid)
            total -= v
    solids = []
    for outer, holes, _, _ in bodies:
        if not holes:
            solids.append(outer)
            continue
        maker = BRepBuilderAPI_MakeSolid(TopoDS.Shell(TopExp_Explorer(outer, TopAbs_SHELL).Current()))
        for h in holes:
            inner = TopoDS.Shell(TopExp_Explorer(h, TopAbs_SHELL).Current())
            maker.Add(TopoDS.Shell(inner.Reversed()))
        solids.append(fixed(maker.Solid()))
    # (bodies side by side: fusing the ones that overlap is left until the part has
    # checked out, see fuse_overlapping; on an assembly of many bodies it is slow)
    solid = solids[0] if len(solids) == 1 else compound(solids)
    return solid, len(solids), len(cavities), total


def inward_shells(mesh):
    """Does any closed shell of the mesh face inward (a cavity)? Without one, no shell
    built from it can be a cavity, and sorting cavities from bodies (slow when a part is
    many bodies) is skipped. A shell here: triangles joined through edges they share."""
    if "inward" not in mesh.__dict__:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        T = mesh.tris
        e = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
        owner = np.tile(np.arange(len(T)), 3)
        key = e[:, 0].astype(np.int64) * (len(mesh.pts) + 1) + e[:, 1]
        order = np.argsort(key, kind="stable")
        k = key[order]
        same = np.nonzero(k[1:] == k[:-1])[0]
        a, b = owner[order[same]], owner[order[same + 1]]
        graph = coo_matrix((np.ones(len(a)), (a, b)), shape=(len(T), len(T)))
        n, label = connected_components(graph, directed=False)
        p, q, r = (mesh.pts[T[:, i]] for i in range(3))
        vol = np.bincount(label, np.einsum("ij,ij->i", p, np.cross(q, r)), n)
        mesh.inward = bool((vol < 0).any())
    return mesh.inward


def fuse_overlapping(shape):
    """Fuse a checked part's bodies where they overlap (a rib sunk into a floor). Kept
    only if the result is a valid solid whose volume lies between the biggest body's
    and all of them together; otherwise they are handed over side by side."""
    solids, ex = [], TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        solids.append(ex.Current())
        ex.Next()
    if len(solids) < 2:
        return shape
    vols = [volume(s) for s in solids]
    slack = UNIFY_VOLUME * sum(vols) + 1e-3
    fused = fuse_checked(solids[:1], solids[1:], max(vols) - slack, sum(vols) + slack)
    return shape if fused is None else fused


def fuse_checked(args, tools, low, high):
    """The fuse of the shapes if it is a valid solid with a volume between low and high,
    else None. Done in a worker process when there are any, and given up after
    FUSE_SECONDS (a boolean on a part of many thousand faces can grind on for minutes;
    bodies merely touching or apart gain nothing from it anyway)."""
    pool = workers.get()
    if pool is None:
        return _fuse_checked(args, tools, low, high)
    import shutil
    import tempfile
    from concurrent.futures import TimeoutError as Timeout
    tmp = tempfile.mkdtemp(prefix="stl2curves_")
    try:
        a, t, out = (os.path.join(tmp, n) for n in ("args.brep", "tools.brep", "fused.brep"))
        BRepTools.Write_s(compound(args), a)
        BRepTools.Write_s(compound(tools), t)
        job = pool.submit(_fuse_job, (a, t, low, high, out))
        try:
            if not job.result(timeout=FUSE_SECONDS):
                return None
        except Timeout:
            workers.abandoned()     # (the worker is left to finish; then restarted)
            return None
        except Exception:
            workers.broken()
            return None
        fused = TopoDS_Shape()
        BRepTools.Read_s(fused, out, BRep_Builder())
        return fused
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _fuse_job(job):
    """In a worker process: _fuse_checked on shapes read from BRep files, the result
    written to one. True if there is a result."""
    a, t, low, high, out = job
    shapes = []
    for path in (a, t):
        comp = TopoDS_Shape()
        BRepTools.Read_s(comp, path, BRep_Builder())
        it, items = TopoDS_Iterator(comp), []
        while it.More():
            items.append(it.Value())
            it.Next()
        shapes.append(items)
    fused = _fuse_checked(shapes[0], shapes[1], low, high)
    return fused is not None and BRepTools.Write_s(fused, out)


def _fuse_checked(args, tools, low, high):
    try:
        fused = boolean(BRepAlgoAPI_Fuse, args, tools)
    except RuntimeError:
        return None
    if count(fused, TopAbs_FACE) and BRepCheck_Analyzer(fused).IsValid() and low <= volume(fused) <= high:
        return fused
    return None


def attempt(mesh, features, mesh_tol, tol, fuse, faceted_volume):
    """Build and sew the part with these features. Returns (result, to_drop): result is
    (shape, bodies, cavities), or None if it fails the checks; to_drop lists the patches
    to leave out next time (faces that couldn't be built, or those blamed for a gap or an
    invalid face), and attempt.blamed says whether they were only blamed. attempt.why
    names the check that failed ("build", "free", "nosolid", "bodies", "volume",
    "volume2", "invalid"; "ok"), for diagnosing a slow or poor conversion."""
    attempt.blamed = False
    attempt.why = "ok"
    comp, shells, failed = build_faces(mesh, features, mesh_tol)
    if failed:
        attempt.why = "build"
        return None, failed
    attempt.blamed = True
    worst = max((f.worst for f in features), default=0.0)
    sewn, free = sew(comp, max(tol, min(0.2, 1.5 * worst)), shells)
    if free:
        # patches whose faces left gaps: drop just those and try again
        attempt.why = "free"
        return None, culprits(mesh, features, free)
    try:
        shape, nb, nv, signed = solids_from_shells(sewn, inward_shells(mesh))
    except RuntimeError:
        attempt.why = "nosolid"
        return None, []
    _mesh_defective(mesh, tol, fuse)
    if mesh.bare_bodies and nb > mesh.bare_bodies:
        attempt.why = "bodies"
        return None, []     # a body the bare facets don't make: a face ran off somewhere
    change = sum(f.change for f in features)
    expected = faceted_volume + change
    allowed = sum(f.tolerance for f in features) + 1e-6 * abs(faceted_volume) + 1e-3
    if abs(signed - expected) > allowed:
        # (OCC's integration can be off where the sewing closed wide gaps: measured
        # again from the faces' tessellation, which follows the edges they share)
        tess = tessellated_volume(shape)
        if tess is not None and abs(tess - expected) <= allowed:
            signed = tess
    if abs(signed - expected) > allowed:
        # a mesh whose bare facets already don't add up to its volume (it crosses itself,
        # or is many bodies touching at corners, some too thin to make a solid) is
        # measured against what the bare facets give instead
        drift = abs(mesh.bare_volume - faceted_volume) if mesh.bare_volume is not None else None
        defective = _mesh_defective(mesh, tol, fuse)
        # (a sound mesh only by what flattening its facets moves: anything more would be
        # the bare build itself going wrong, and is no reference)
        if drift is None or (not defective and drift > BARE_DRIFT * abs(faceted_volume) + 1e-3):
            attempt.why = "volume"
            return None, []
        expected = mesh.bare_volume + change
        # (a defective mesh give or take what the defect itself does: it sews a little
        # differently each time)
        slack = drift if defective else 0.0
        if abs(signed - expected) > allowed + slack:
            attempt.why = "volume2"
            return None, []
    check = BRepCheck_Analyzer(shape)
    if not check.IsValid():
        fixed_shape = None
        if not mesh.__dict__.get("defective"):      # (no repair mends a mesh that crosses itself)
            # (a straight edge running on tangent into an arc can pass for a crossing at
            # the default precision; a coarser one, still far below the sewing, mends it)
            for precision in (None,) + FIX_PRECISIONS:
                fix = ShapeFix_Shape(shape)
                if precision:
                    fix.SetPrecision(precision)
                fix.Perform()
                fixed_shape = fix.Shape()
                if BRepCheck_Analyzer(fixed_shape).IsValid() and abs(volume(fixed_shape) - expected) <= allowed:
                    break
        if (fixed_shape is not None and BRepCheck_Analyzer(fixed_shape).IsValid()
                and abs(volume(fixed_shape) - expected) <= allowed):
            shape = fixed_shape
        else:
            # patches next to faces that came out invalid: drop just those and try again
            blame = culprits(mesh, features, invalid_face_points(shape, check))
            if blame or not _mesh_defective(mesh, tol, fuse):
                attempt.why = "invalid"
                return None, blame
            # nothing to blame and the mesh as bare facets fails the same way: the
            # defect is in the mesh itself (it touches itself, say), not in the curves
    trial = mesh.__dict__.get("file_trial", False)
    if (features and not mesh.__dict__.get("file_off")
            and (trial or mesh.__dict__.get("file_rounds", 0) < FILE_ROUNDS)):
        # The solid is only as good as its STEP file. A file keeps no tolerances, and a
        # reader works them out again from the geometry, tight: a face that checked out
        # only within the wide tolerance sewing gave its edges can come back crossing
        # itself, split up, or (an edge closed up by sewing) running round a whole circle.
        trouble = file_trouble(shape, allowed)
        if trouble:
            blame = culprits(mesh, features, trouble)
            if blame:
                # (a few rounds at most: trouble the patches nearby don't cure is left to
                # the final check's warning rather than costing every patch round it)
                if not trial:
                    mesh.file_rounds = mesh.__dict__.get("file_rounds", 0) + 1
                attempt.why = "file"
                return None, blame
            # (nothing near to blame: the bare facets themselves, kept as they are)
    return (shape, nb, nv), []


def read_step(path):
    """The shape in a STEP file (OpenCascade's transfer statistics silenced)."""
    sys.stdout.flush()
    saved, devnull = os.dup(1), os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    try:
        r = STEPControl_Reader()
        r.ReadFile(str(path))
        r.TransferRoots()
        return r.OneShape()
    finally:
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


def file_trouble(shape, allowed):
    """Write the shape to STEP and read it back: None if it comes back a valid solid of
    the same volume (within allowed, mm^3), else points where it doesn't (maybe none
    found). Only real changes count: a face that comes back invalid, goes missing or
    changes its area by more than FILE_AREA (small faces read back a percent or so off,
    their tolerant edges laid out afresh, and that is harmless); and on a big face only
    the outline (wire) that changed, not every hole in it."""
    tmp = tempfile.mkdtemp(prefix="stl2curves_")
    try:
        path = os.path.join(tmp, "check.step")
        write_step(shape, path)
        back = read_step(path)
    except Exception:
        return []
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    check = BRepCheck_Analyzer(back)
    # (a reader may split a face in two: harmless, so faces aren't counted)
    if check.IsValid() and abs(volume(back) - volume(shape)) <= allowed:
        return None

    def faces(s):
        out, ex = [], TopExp_Explorer(s, TopAbs_FACE)
        while ex.More():
            face = TopoDS.Face(ex.Current())
            ex.Next()
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, props)
            c = props.CentreOfMass()
            out.append((face, int(BRepAdaptor_Surface(face).GetType()), props.Mass(), np.array([c.X(), c.Y(), c.Z()])))
        return out

    def partner(item, pool, table):
        # the face in the other shape most like this one: same kind of surface, middle
        # nearby (a big face's middle moves further when it changes), area closest
        _, kind, area, centre = item
        kinds, areas, centres = table
        reach = FILE_MATCH + 0.02 * math.sqrt(abs(area))
        near = np.nonzero((kinds == kind) & (np.linalg.norm(centres - centre, axis=1) <= reach))[0]
        return pool[near[np.argmin(np.abs(areas[near] - area))]] if len(near) else None

    def table(pool):
        return (np.array([x[1] for x in pool]), np.array([x[2] for x in pool]),
                np.array([x[3] for x in pool]).reshape(-1, 3))

    mine, theirs = faces(shape), faces(back)
    mine_t, theirs_t = table(mine), table(theirs)
    points = []
    for item in mine:
        other = partner(item, theirs, theirs_t)
        if other is None or abs(other[2] - item[2]) > FILE_AREA + 0.005 * abs(item[2]):
            points += _changed_outline(item[0], None if other is None else other[0])
    for item in theirs:
        if not check.IsValid(item[0]):
            other = partner(item, mine, mine_t)
            points += (_changed_outline(other[0], item[0]) if other is not None
                       else _vertex_points(item[0]))
    return points or invalid_face_points(back, check)


def _vertex_points(shape):
    out, vx = [], TopExp_Explorer(shape, TopAbs_VERTEX)
    while vx.More():
        p = BRep_Tool.Pnt_s(TopoDS.Vertex(vx.Current()))
        out.append((p.X(), p.Y(), p.Z()))
        vx.Next()
    return out


def _outlines(face):
    """(wire, length, middle, vector area) of each of the face's wires, from points along
    its edges in order (the vector area turns round when the wire does)."""
    out, wx = [], TopExp_Explorer(face, TopAbs_WIRE)
    while wx.More():
        wire = TopoDS.Wire(wx.Current())
        wx.Next()
        pts, we = [], BRepTools_WireExplorer(wire, face)
        while we.More():
            edge = we.Current()
            we.Next()
            c = BRepAdaptor_Curve(edge)
            t = np.linspace(c.FirstParameter(), c.LastParameter(), 9)
            if edge.Orientation() == TopAbs_REVERSED:
                t = t[::-1]
            pts += [(p.X(), p.Y(), p.Z()) for p in (c.Value(x) for x in t[:-1])]
        if not pts:
            continue
        P = np.array(pts)
        Q = np.roll(P, -1, axis=0)
        step = np.linalg.norm(Q - P, axis=1)
        # (the middle weighted by length: a reader may split an edge in two)
        middle = ((P + Q) / 2 * step[:, None]).sum(axis=0) / max(step.sum(), 1e-12)
        out.append((wire, float(step.sum()), middle, 0.5 * np.cross(P, Q).sum(axis=0)))
    return out


def _changed_outline(face, other):
    """Points on each wire of the face that the other face (its version read back from
    the file, or None) doesn't have alike; all the face's corners if every wire matches."""
    theirs = _outlines(other) if other is not None else []
    out = []
    for wire, length, middle, turn in _outlines(face):
        if not any(abs(l - length) <= FILE_MATCH * (1 + length) and np.linalg.norm(m - middle) <= FILE_MATCH
                   and np.linalg.norm(t - turn) <= FILE_MATCH * (1 + np.linalg.norm(turn))
                   for _, l, m, t in theirs):
            out += _vertex_points(wire)
    return out or _vertex_points(face)


def _mesh_defective(mesh, tol, fuse):
    """Is the mesh, built from bare facets, already not a valid solid? Sets
    mesh.bare_volume and mesh.bare_bodies to the volume and number of bodies the bare
    facets give (None and 0 if they make no solid)."""
    if "defective" not in mesh.__dict__:
        got, ahead = None, mesh.__dict__.pop("bare_ahead", None)
        if ahead is not None:
            try:
                got = ahead.result()
            except Exception:
                got = None
        mesh.defective, mesh.bare_volume, mesh.bare_bodies = got or _bare(mesh, tol)
    return mesh.defective


def _bare(mesh, tol):
    """(defective, volume, bodies) of the mesh built from bare facets."""
    comp, shells, _ = build_faces(mesh, [], TOL)
    sewn, free = sew(comp, tol, shells)
    try:
        shape, bodies, _, vol = solids_from_shells(sewn, inward_shells(mesh))
        return bool(free) or not BRepCheck_Analyzer(shape).IsValid(), vol, bodies
    except RuntimeError:
        return True, None, 0


def _bare_job(job):
    """In a worker process: _bare on the shared mesh."""
    key, tol = job
    return _bare(workers.load(key), tol)


def _bare_ahead(mesh, tol):
    """Start building a big mesh from bare facets in a worker process (_mesh_defective
    picks it up), so it is ready by the time the first attempt is checked."""
    pool = workers.get() if len(mesh.farea) >= BARE_AHEAD else None
    if pool is not None and "defective" not in mesh.__dict__:
        try:
            if "shared_as" not in mesh.__dict__:
                mesh.shared_as = workers.share(mesh)
            mesh.bare_ahead = pool.submit(_bare_job, (mesh.shared_as, tol))
        except Exception:
            mesh.__dict__.pop("bare_ahead", None)


def culprits(mesh, features, points):
    """Which patches to drop for trouble at these points: the smooth blends right there if
    any (they give back the exact pieces they replaced), else just the patch nearest to
    each trouble spot (the next attempt shows whether that was enough)."""
    near = features_near(mesh, features, points)
    if not near:
        return []
    P = np.asarray(points)
    dist = {}
    for k in near:
        dist[k] = point_facet_distance(mesh, features[k].facets, P)     # per trouble spot
    # blends are the likeliest cause and the cheapest loss: blame one nearby first
    blends = [k for k in near if features[k].kind == "blend"]
    out = set()
    for i in range(len(P)):
        pool = [k for k in blends if dist[k][i] <= 0.5] or near
        k = min(pool, key=lambda k: dist[k][i])
        if dist[k][i] <= 0.5:
            out.add(k)
    return sorted(out)


def invalid_face_points(shape, check=None):
    """A point on each face of the shape that fails the validity check (check: the
    shape's BRepCheck_Analyzer, if already run: it has every face's verdict)."""
    out, ex = [], TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        face = ex.Current()
        bad = not (check.IsValid(face) if check is not None else BRepCheck_Analyzer(face).IsValid())
        if bad:
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, props)
            c = props.CentreOfMass()
            out.append((c.X(), c.Y(), c.Z()))
            vx = TopExp_Explorer(face, TopAbs_VERTEX)
            while vx.More():
                p = BRep_Tool.Pnt_s(TopoDS.Vertex(vx.Current()))
                out.append((p.X(), p.Y(), p.Z()))
                vx.Next()
        ex.Next()
    return out + misjoined_edge_points(shape)


def misjoined_edge_points(shape):
    """The middle of each edge its shell uses wrongly: run the same way by both faces
    beside it (one of them is flipped), or shared by more than two faces. Every face
    can check out on its own while its shell doesn't; this says where."""
    out, ex = [], TopExp_Explorer(shape, TopAbs_SHELL)
    while ex.More():
        edges = IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher()
        TopExp.MapShapes_s(ex.Current(), TopAbs_EDGE, edges)
        uses = np.zeros((edges.Extent() + 1, 2), int)       # forward, reversed
        fx = TopExp_Explorer(ex.Current(), TopAbs_FACE)
        while fx.More():
            ee = TopExp_Explorer(fx.Current(), TopAbs_EDGE)
            while ee.More():
                o = ee.Current().Orientation()
                if o in (TopAbs_FORWARD, TopAbs_REVERSED):
                    uses[edges.FindIndex(ee.Current()), int(o == TopAbs_REVERSED)] += 1
                ee.Next()
            fx.Next()
        n = uses.sum(axis=1)
        for i in np.nonzero((n > 2) | ((n == 2) & (uses[:, 0] != 1)))[0]:
            edge = TopoDS.Edge(edges.FindKey(int(i)))
            if not BRep_Tool.Degenerated_s(edge):
                c = BRepAdaptor_Curve(edge)
                p = c.Value((c.FirstParameter() + c.LastParameter()) / 2)
                out.append((p.X(), p.Y(), p.Z()))
        ex.Next()
    return out


def _build(mesh, features, mesh_tol, tol, fuse, info):
    """Build and check one group of bodies, leaving faceted any feature that spoils it.
    Returns (shape, bodies, cavities)."""
    faceted_volume = mesh_volume(mesh)
    _bare_ahead(mesh, tol)
    features = settle_blends(mesh, features, mesh_tol, split_blend, info["skipped"])

    def good(subset):
        result, failed = attempt(mesh, subset, mesh_tol, tol, fuse, faceted_volume)
        blamed = []
        while failed:  # patches whose face couldn't be built: drop them and retry
            # a smooth blend whose face won't fit is cut in two and tried again; failing
            # that (or if it was only blamed), it gives back the pieces it replaced
            back = []
            by_file = attempt.why == "file"
            for k in failed:
                f = subset[k]
                rings = half_rings(mesh, f)
                if rings:
                    # a whole ring cut to shape often fails where its outline crosses
                    # the surface's seam; the same surface as two half rings doesn't
                    back += rings
                    continue
                if attempt.blamed:
                    halves = None
                elif f.kind == "blend":
                    halves = split_blend(mesh, f)
                else:
                    # an exact surface whose face couldn't be built: try a smooth blend
                    # over its facets rather than leaving them flat
                    halves = [as_blend(mesh, f.facets, ())] if len(f.facets) >= 3 else None
                if halves:
                    back += halves
                    if f.kind != "blend":
                        info["skipped"].append(f)
                else:
                    info["skipped"].append(f)
                    back += list(f.parts)
                    if by_file:
                        mesh.__dict__.setdefault("file_dropped", []).append(f)
                    if attempt.blamed:
                        blamed.append(f)
            subset = [f for k, f in enumerate(subset) if k not in failed] + back
            result, failed = attempt(mesh, subset, mesh_tol, tol, fuse, faceted_volume)

        # Patches dropped for trouble nearby may have been innocent: once the part
        # builds, give them a second chance, all together first, halving on failure.
        def again(group):
            nonlocal result, subset
            if not group or time_left() < 0:
                return
            trial = [f for f in subset if all(f is not p for b in group for p in b.parts)] + group
            # (the STEP file checked whatever rounds are left: a patch dropped for it must
            # not come back unchecked)
            mesh.file_trial = True
            try:
                got, drop = attempt(mesh, trial, mesh_tol, tol, fuse, faceted_volume)
            finally:
                mesh.file_trial = False
            if got is not None and not drop:
                result, subset = got, trial
                info["skipped"] = [f for f in info["skipped"] if all(f is not b for b in group)]
            elif len(group) > 1:
                again(group[:len(group) // 2])
                again(group[len(group) // 2:])

        if result is not None:
            again(blamed)
        return result, subset

    result, used = good(features)
    if result is None and used:
        # Find the troublemakers by halving: keep every half that builds cleanly.
        accepted = []
        start = time.time()

        def search(group):
            nonlocal result
            if not group:
                return
            if time.time() - start > SEARCH_SECONDS or time_left() < 0:
                info["skipped"] += group      # out of time: leave the rest faceted
                return
            trial, kept = good(accepted + group)
            if trial is not None:
                accepted[:] = kept
                result = trial
                return
            if len(group) == 1:
                info["skipped"].append(group[0])
                return
            half = len(group) // 2
            search(group[:half])
            search(group[half:])

        search(used)
        used = accepted
        if result is None:
            result, used = good([])
    if result is None:
        raise RuntimeError("could not build a closed solid from this mesh")
    dropped = [f for f in mesh.__dict__.pop("file_dropped", []) if all(f is not u for u in used)]
    if dropped:
        allowed = sum(f.tolerance for f in used) + 1e-6 * abs(faceted_volume) + 1e-3
        if file_trouble(result[0], allowed) is not None:
            # dropping the patches blamed for the STEP file didn't mend it: better have
            # them back (the final check warns about the file either way)
            mesh.file_off = True
            skipped = list(info["skipped"])
            info["skipped"] = [f for f in skipped if all(f is not b for b in dropped)]
            got, kept = good([f for f in used if all(f is not p for b in dropped for p in b.parts)] + dropped)
            if got is not None:
                result, used = got, kept
            else:
                info["skipped"] = skipped
    info["restored"] += used
    shape, nb, nv = result
    if fuse and nb > 1:
        shape = fuse_overlapping(shape)
        nb = count(shape, TopAbs_SOLID)
    return shape, nb, nv


def stl_to_solid(*args, time_limit=TIME_LIMIT, **kwargs):
    """Convert one STL: (shape, info). See _stl_to_solid. time_limit (seconds, None for
    none): fitting blends and hunting down troublemakers stop at their share of it, and
    whatever has checked out by then is handed over."""
    features_mod._budget = time_limit
    features_mod._deadline = time.time() + time_limit if time_limit else None

    try:
        return _stl_to_solid(*args, **kwargs)
    finally:
        features_mod._deadline = features_mod._budget = None
        workers.finish()


def _stl_to_solid(path, tol, fuse=True, curves=True, true_size=False, blends=True, mend=True,
                  simplify_to=None):
    """simplify_to: how far (mm) thinning out an over-dense mesh may move its surface;
    None: only above AUTO_SIMPLIFY triangles, at SIMPLIFY_ERROR; 0: never."""
    # (a mesh file, or (points, triangles) or (points, triangles, part of each triangle):
    # a 3MF object's parts can overlap, and are each mended and built on their own, then
    # joined by a boolean union)
    part = None
    if isinstance(path, (str, os.PathLike)):
        pts, tris = load_stl(path)
    else:
        pts, tris = path[0], path[1]
        if len(path) > 2 and len(np.unique(path[2])) > 1:
            part = np.asarray(path[2])
    if len(tris) >= WORKERS_FROM:
        workers.start(len(tris))    # (worker processes, started while the mesh is repaired)
    info = {"triangles": len(tris), "restored": [], "skipped": [], "size": None, "snapped": 0,
            "repairs": [], "simplified": None}
    if mend:
        if part is None:
            pts, tris, info["repairs"] = repair(pts, tris)
        else:
            pts, tris, part, info["repairs"] = _repair_parts(pts, tris, part)
    if simplify_to is None:
        simplify_to = SIMPLIFY_ERROR if len(tris) > AUTO_SIMPLIFY else 0
    if simplify_to:
        before = len(tris)
        if part is None:
            tris = simplify(pts, tris, simplify_to)
        else:
            labels = np.unique(part)
            pieces = [simplify(pts, tris[part == p], simplify_to) for p in labels]
            part = np.concatenate([np.full(len(t), p) for p, t in zip(labels, pieces)])
            tris = np.vstack(pieces)
        info["simplified"] = (before, len(tris), simplify_to)
    # bodies touching at an edge are built apart (and so are a 3MF object's parts)
    groups = (split_bodies(pts, tris) if part is None
              else [g for p in np.unique(part) for g in split_bodies(pts, tris[part == p])])
    if curves:
        parts = [analyze(pts, g) for g in groups]
        info["size"] = guess = guess_size([p[:2] for p in parts])
        round_unit = None
        if true_size and guess is not None:
            round_unit = guess.unit
            if guess.factor != 1.0:
                parts = [analyze(pts * guess.factor, g) for g in groups]
        for mesh, features, _ in parts:
            info["snapped"] += snap(mesh, features, round_unit)
        if blends:
            parts = [(mesh, add_blends(mesh, features), t) for mesh, features, t in parts]
    else:
        parts = [(Mesh(pts, g), [], TOL) for g in groups]
    shapes, nb, nv = [], 0, 0
    for mesh, features, mesh_tol in parts:
        shape, b, v = _build(mesh, features, mesh_tol, tol, fuse, info)
        shapes.append(shape)
        nb, nv = nb + b, nv + v
        info["mesh_defects"] = info.get("mesh_defects", False) or mesh.__dict__.get("defective", False)
    shape = shapes[0]
    if len(shapes) > 1:
        # bodies that met at an edge were built apart; join them now, unless the join
        # comes out invalid or loses material (the boolean can return nothing at all
        # where bodies only touch): then they are handed over side by side
        shape = compound(shapes)
        if fuse:
            apart = sum(volume(s) for s in shapes)
            slack = UNIFY_VOLUME * abs(apart) + 1e-3
            # (a 3MF object's parts may overlap: then the union is less than their sum,
            # but no less than the biggest)
            low = apart if part is None else max(volume(s) for s in shapes)
            fused = fuse_checked(shapes[:1], shapes[1:], low - slack, apart + slack)
            if fused is not None:
                shape = fused
        nb = count(shape, TopAbs_SOLID)

    # Tidy up, but only keep the result if it is still a valid solid of the same volume;
    # otherwise hand over the checked shape as it was built. First faces on one surface
    # are merged into one (a flat face whose outline defeated the face builder, a step or
    # ridge a thousandth of a millimetre high in an export, was built triangle by
    # triangle), failing that just edges split along one line.
    checked = shape
    target = volume(checked)
    shape = untidy = None
    for faces in (True, False):
        try:
            unify = ShapeUpgrade_UnifySameDomain(checked, True, faces, False)
            unify.SetLinearTolerance(UNIFY_TOL)
            unify.SetAngularTolerance(math.radians(UNIFY_DEG))
            unify.Build()
            sf = ShapeFix_Shape(unify.Shape())
            sf.Perform()
            tidy = sf.Shape()
        except Exception:
            continue
        if abs(volume(tidy) - target) <= UNIFY_VOLUME * abs(target) + 1e-3:
            if BRepCheck_Analyzer(tidy).IsValid():
                shape = tidy
                break
            untidy = untidy or tidy
    if shape is None:
        shape = checked
    if not BRepCheck_Analyzer(shape).IsValid():
        if BRepCheck_Analyzer(checked).IsValid():
            shape = checked
        else:
            sf = ShapeFix_Shape(checked)
            sf.Perform()
            if BRepCheck_Analyzer(sf.Shape()).IsValid():
                shape = sf.Shape()
            elif untidy is not None:
                shape = untidy      # (no worse than the shape as checked, and tidier)
    # (and the STEP file must read back as built: the tidied shape failing that, the
    # shape as checked)
    slack = FILE_VOLUME * abs(target) + 1e-3
    trouble = file_trouble(shape, slack)
    if trouble is not None and shape is not checked and file_trouble(checked, slack) is None:
        shape, trouble = checked, None
    info["file_ok"] = trouble is None

    info["bodies"], info["voids"] = nb, nv
    info["faces"] = count(shape, TopAbs_FACE)
    if not info["faces"]:
        raise RuntimeError("the solid came out empty")
    info["volume"] = volume(shape)
    info["valid"] = BRepCheck_Analyzer(shape).IsValid()
    return shape, info


def _repair_parts(pts, tris, part):
    """repair() on each part of a 3MF object on its own (two parts pressed together are
    no wall of zero thickness to cut out, nor overlapping ones a mesh crossing itself):
    (points, triangles, part of each triangle, what was done)."""
    P, T, L, done = [], [], [], {}
    base = 0
    for p in np.unique(part):
        used, inv = np.unique(tris[part == p], return_inverse=True)
        q, t, notes = repair(pts[used], inv.reshape(-1, 3))
        P.append(q)
        T.append(t + base)
        L.append(np.full(len(t), p))
        base += len(q)
        for note in notes:
            # (the same repair in two parts is reported once, with the counts added up)
            head, _, rest = note.partition(" ")
            if head.isdigit():
                done[rest] = done.get(rest, 0) + int(head)
            else:
                done[note] = None
    notes = [note if n is None else f"{n} {note}" for note, n in done.items()]
    return np.vstack(P), np.vstack(T), np.concatenate(L), notes


def write_step(shape, path):
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    # OpenCascade prints a wall of transfer statistics to stdout; silence it.
    sys.stdout.flush()
    saved, devnull = os.dup(1), os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    try:
        w = STEPControl_Writer()
        w.Transfer(shape, STEPControl_AsIs)
        status = w.Write(str(path))
    finally:
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)
    if status != IFSelect_RetDone:
        raise RuntimeError(f"failed writing {path}")


def describe_features(info, details=False):
    text = ""
    if info["restored"]:
        text = f"  rebuilt as true curves: {summarize(info['restored'])}"
    if info["skipped"]:
        text += f"\n  left faceted (failed checks): {summarize(info['skipped'])}"
    if details:
        for f in info["restored"]:
            text += f"\n    {f.describe()}"
        for f in info["skipped"]:
            text += f"\n    SKIPPED {f.describe()}"
    return text.lstrip("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--out", help="output folder (default: next to each STL)")
    ap.add_argument("--merge", help="also write all parts into one STEP with this name")
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--no-fuse", action="store_true",
                    help="keep separate/overlapping bodies in a file separate instead of unioning them")
    ap.add_argument("--no-curves", action="store_true",
                    help="don't rebuild curved areas as true curves")
    ap.add_argument("--details", action="store_true", help="list every rebuilt feature")
    ap.add_argument("--no-blends", action="store_true",
                    help="don't turn curved areas no simple surface fits into smooth freeform faces")
    ap.add_argument("--no-repair", action="store_true",
                    help="don't mend mesh defects (duplicates, slivers, cracks, holes, flipped or crossing triangles)")
    ap.add_argument("--simplify", type=float, metavar="MM",
                    help=f"thin out the mesh first, moving its surface by at most MM (0: never; default: "
                         f"{SIMPLIFY_ERROR} mm for meshes over {AUTO_SIMPLIFY:,} triangles)")
    ap.add_argument("--time-limit", type=float, metavar="SECONDS", default=TIME_LIMIT,
                    help=f"stop fitting blends and hunting down troublemakers after about this long "
                         f"and keep what checks out (0: no limit; default {TIME_LIMIT})")
    ap.add_argument("--true-size", action="store_true",
                    help="rebuild at the apparent design size, with radii snapped to round values")
    args = ap.parse_args()

    files = []
    for p in map(Path, args.inputs):
        files += sorted(f for f in p.iterdir() if f.suffix.lower() in (".stl", ".3mf")) if p.is_dir() else [p]
    if not files:
        sys.exit("No STL or 3MF files found.")
    # (a 3MF file holds any number of objects: each is converted to a STEP of its own)
    jobs = []
    for f in files:
        if f.suffix.lower() == ".3mf":
            objects = read_3mf(f)
            if not objects:
                print(f"{f.name}: no objects to convert", flush=True)
            jobs += [(f"{f.name}: {name}", (pts, tris, part), f, f"{f.stem} - {name}")
                     for name, pts, tris, part in objects]
        else:
            jobs.append((f.name, f, f, f.stem))

    shapes = []
    for label, source, f, stem in jobs:
        t = time.time()
        print(f"{label}: converting...", flush=True)
        try:
            shape, info = stl_to_solid(source, args.tol, not args.no_fuse, not args.no_curves, args.true_size,
                                       not args.no_blends, not args.no_repair, args.simplify,
                                       time_limit=args.time_limit or None)
        except Exception as e:
            # (one part that can't be converted shouldn't stop the rest of the batch)
            print(f"  FAILED: {e} ({time.time() - t:.1f}s); no STEP written", flush=True)
            continue
        out_dir = Path(args.out) if args.out else f.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / (stem + ".step")
        write_step(shape, out)
        shapes.append(shape)
        nb, nv = info["bodies"], info["voids"]
        extra = f" ({nb} bodies{', ' + str(nv) + ' cavities' if nv else ''})" if nb + nv > 1 else ""
        print(f"  {info['triangles']} triangles -> {info['faces']} faces{extra}, "
              f"volume {info['volume']:,.1f} mm^3, "
              f"{'valid solid' if info['valid'] else 'WARNING: check geometry' + (' (the STL itself is not a clean solid: it touches or crosses itself)' if info.get('mesh_defects') else '')}, "
              f"{'' if info.get('file_ok', True) else 'WARNING: the STEP file does not read back exactly as built, '}"
              f"{time.time() - t:.1f}s -> {out}")
        if info["repairs"]:
            print("  mended the mesh: " + "; ".join(info["repairs"]))
        if info["simplified"]:
            before, after, err = info["simplified"]
            print(f"  thinned out: {before:,} -> {after:,} triangles (surface moved {err} mm at most)")
        if info["restored"] or info["skipped"]:
            print(describe_features(info, args.details))
        if info["size"] is not None:
            note = info["size"].describe()
            if args.true_size and info["size"].factor != 1.0:
                note = note.split(" Use --true-size")[0] + " Rebuilt at that size."
            print("  " + note)

    if args.merge:
        comp = compound(shapes)
        out_dir = Path(args.out) if args.out else files[0].parent
        out = out_dir / (Path(args.merge).stem + ".step")
        write_step(comp, out)
        print(f"Combined file -> {out}")


if __name__ == "__main__":
    main()
