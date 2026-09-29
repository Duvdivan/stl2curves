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
import os
import sys
import time
from pathlib import Path

import numpy as np
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_IN
from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.IFSelect import IFSelect_RetDone
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.collections import List_TopoDS_Shape

from features import load_stl, analyze, summarize, snap, half_rings, Mesh, TOL
from sizing import guess_size
from build import build_faces, sew, features_near, settle_blends
from bodies import split_bodies
from blends import add_blends, split as split_blend, _blend as as_blend
from repair import repair
from simplify import simplify

SEARCH_SECONDS = 600    # time allowed for hunting down patches that spoil the solid
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


def solids_from_shells(sewn, fuse):
    """Turn sewn shells into solids. A shell inside another is a cavity, the rest are bodies.

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
    shells.sort(key=lambda x: -x[2])        # biggest first: an enclosing shell comes first

    def points_of(shape, most=24):
        # corners and face centres (a rib spanning wall to wall has every corner buried
        # in the walls, but the middle of its top face is out in the open)
        pts, ex = [], TopExp_Explorer(shape, TopAbs_FACE)
        while ex.More():
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(ex.Current(), props)
            pts.append(props.CentreOfMass())
            ex.Next()
        ex = TopExp_Explorer(shape, TopAbs_VERTEX)
        while ex.More():
            pts.append(BRep_Tool.Pnt_s(TopoDS.Vertex(ex.Current())))
            ex.Next()
        return pts[::max(1, len(pts) // most)]

    def inside(solid, pts):
        return all(BRepClass3d_SolidClassifier(solid, p, 1e-6).State() == TopAbs_IN for p in pts)

    # A cavity lies wholly inside its host; an overlapping body (a rib sunk into a
    # floor, say) has some corners outside, and gets fused on instead.
    bodies, cavities, total = [], [], 0.0      # bodies: [outer solid, [cavity solids]]
    for shell, solid, v in shells:
        pts = points_of(shell)
        host = next((b for b in bodies if inside(b[0], pts)), None)
        if host is None:
            bodies.append([solid, []])
            total += v
        else:
            host[1].append(solid)
            cavities.append(solid)
            total -= v
    solids = []
    for outer, holes in bodies:
        if not holes:
            solids.append(outer)
            continue
        maker = BRepBuilderAPI_MakeSolid(TopoDS.Shell(TopExp_Explorer(outer, TopAbs_SHELL).Current()))
        for h in holes:
            inner = TopoDS.Shell(TopExp_Explorer(h, TopAbs_SHELL).Current())
            maker.Add(TopoDS.Shell(inner.Reversed()))
        solids.append(fixed(maker.Solid()))
    solid = solids[0]
    if fuse and len(solids) > 1:
        solid = boolean(BRepAlgoAPI_Fuse, solids[:1], solids[1:])
    elif len(solids) > 1:
        solid = compound(solids)
    return solid, len(solids), len(cavities), total


def attempt(mesh, features, mesh_tol, tol, fuse, faceted_volume):
    """Build and sew the part with these features. Returns (result, to_drop): result is
    (shape, bodies, cavities), or None if it fails the checks; to_drop lists the patches
    to leave out next time (faces that couldn't be built, or those blamed for a gap or an
    invalid face), and attempt.blamed says whether they were only blamed."""
    attempt.blamed = False
    comp, shells, failed = build_faces(mesh, features, mesh_tol)
    if failed:
        return None, failed
    attempt.blamed = True
    worst = max((f.worst for f in features), default=0.0)
    sewn, free = sew(comp, max(tol, min(0.2, 1.5 * worst)), shells)
    if free:
        # patches whose faces left gaps: drop just those and try again
        return None, culprits(mesh, features, free)
    try:
        shape, nb, nv, signed = solids_from_shells(sewn, fuse)
    except RuntimeError:
        return None, []
    change = sum(f.change for f in features)
    expected = faceted_volume + change
    allowed = sum(f.tolerance for f in features) + 1e-6 * abs(faceted_volume) + 1e-3
    if abs(signed - expected) > allowed:
        # a mesh whose bare facets already don't add up to its volume (it crosses itself,
        # say) is measured against what the bare facets give instead
        if not _mesh_defective(mesh, tol, fuse) or mesh.bare_volume is None:
            return None, []
        # (give or take what the defect itself does: it sews a little differently each time)
        expected = mesh.bare_volume + change
        if abs(signed - expected) > allowed + abs(mesh.bare_volume - faceted_volume):
            return None, []
    if not BRepCheck_Analyzer(shape).IsValid():
        fixed_shape = None
        if not mesh.__dict__.get("defective"):      # (no repair mends a mesh that crosses itself)
            fix = ShapeFix_Shape(shape)
            fix.Perform()
            fixed_shape = fix.Shape()
        if (fixed_shape is not None and BRepCheck_Analyzer(fixed_shape).IsValid()
                and abs(volume(fixed_shape) - expected) <= allowed):
            shape = fixed_shape
        else:
            # patches next to faces that came out invalid: drop just those and try again
            blame = culprits(mesh, features, invalid_face_points(shape))
            if blame or not _mesh_defective(mesh, tol, fuse):
                return None, blame
            # nothing to blame and the mesh as bare facets fails the same way: the
            # defect is in the mesh itself (it touches itself, say), not in the curves
    return (shape, nb, nv), []


def _mesh_defective(mesh, tol, fuse):
    """Is the mesh, built from bare facets, already not a valid solid? Sets
    mesh.bare_volume to the volume the bare facets give."""
    if "defective" not in mesh.__dict__:
        comp, shells, _ = build_faces(mesh, [], TOL)
        sewn, free = sew(comp, tol, shells)
        mesh.bare_volume = None
        try:
            shape, _, _, mesh.bare_volume = solids_from_shells(sewn, fuse)
            mesh.defective = bool(free) or not BRepCheck_Analyzer(shape).IsValid()
        except RuntimeError:
            mesh.defective = True
    return mesh.defective


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
        V = mesh.pts[np.unique(np.concatenate([mesh.fverts[x] for x in features[k].facets]))]
        dist[k] = np.linalg.norm(P[:, None] - V[None], axis=2).min(axis=1)   # per trouble spot
    # blends are the likeliest cause and the cheapest loss: blame one nearby first
    blends = [k for k in near if features[k].kind == "blend"]
    out = set()
    for i in range(len(P)):
        pool = [k for k in blends if dist[k][i] <= 0.5] or near
        k = min(pool, key=lambda k: dist[k][i])
        if dist[k][i] <= 0.5:
            out.add(k)
    return sorted(out)


def invalid_face_points(shape):
    """A point on each face of the shape that fails the validity check."""
    out, ex = [], TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        face = ex.Current()
        if not BRepCheck_Analyzer(face).IsValid():
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
    return out


def _build(mesh, features, mesh_tol, tol, fuse, info):
    """Build and check one group of bodies, leaving faceted any feature that spoils it.
    Returns (shape, bodies, cavities)."""
    faceted_volume = mesh_volume(mesh)
    features = settle_blends(mesh, features, mesh_tol, split_blend, info["skipped"])

    def good(subset):
        result, failed = attempt(mesh, subset, mesh_tol, tol, fuse, faceted_volume)
        blamed = []
        while failed:  # patches whose face couldn't be built: drop them and retry
            if attempt.blamed:
                blamed += [subset[k] for k in failed if subset[k].kind == "blend"]
            # a smooth blend whose face won't fit is cut in two and tried again; failing
            # that (or if it was only blamed), it gives back the pieces it replaced
            back = []
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
            subset = [f for k, f in enumerate(subset) if k not in failed] + back
            result, failed = attempt(mesh, subset, mesh_tol, tol, fuse, faceted_volume)
        # A blend dropped for trouble nearby may have been innocent: once the part
        # builds, give each one a second chance on its own.
        for b in blamed[:6] if result is not None else []:
            trial = [f for f in subset if all(f is not p for p in b.parts)] + [b]
            got, again = attempt(mesh, trial, mesh_tol, tol, fuse, faceted_volume)
            if got is not None and not again:
                result, subset = got, trial
                info["skipped"] = [f for f in info["skipped"] if f is not b]
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
            if time.time() - start > SEARCH_SECONDS:
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
    info["restored"] += used
    return result


def stl_to_solid(path, tol, fuse=True, curves=True, true_size=False, blends=True, mend=True,
                 simplify_to=None):
    """simplify_to: how far (mm) thinning out an over-dense mesh may move its surface;
    None: only above AUTO_SIMPLIFY triangles, at SIMPLIFY_ERROR; 0: never."""
    pts, tris = load_stl(path)
    info = {"triangles": len(tris), "restored": [], "skipped": [], "size": None, "snapped": 0,
            "repairs": [], "simplified": None}
    if mend:
        pts, tris, info["repairs"] = repair(pts, tris)
    if simplify_to is None:
        simplify_to = SIMPLIFY_ERROR if len(tris) > AUTO_SIMPLIFY else 0
    if simplify_to:
        before = len(tris)
        tris = simplify(pts, tris, simplify_to)
        info["simplified"] = (before, len(tris), simplify_to)
    groups = split_bodies(pts, tris)          # bodies touching at an edge are built apart
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
        # bodies that met at an edge were built apart; join them now
        shape = boolean(BRepAlgoAPI_Fuse, shapes[:1], shapes[1:]) if fuse else compound(shapes)
        nb = count(shape, TopAbs_SOLID)

    # Tidy up (merge edges split along one line), but only keep the result if it is still
    # a valid solid; otherwise hand over the checked shape as it was built
    checked = shape
    try:
        unify = ShapeUpgrade_UnifySameDomain(shape, True, False, False)  # merge edges only
        unify.Build()
        sf = ShapeFix_Shape(unify.Shape())
        sf.Perform()
        shape = sf.Shape()
    except Exception:
        shape = checked
    if not BRepCheck_Analyzer(shape).IsValid():
        if BRepCheck_Analyzer(checked).IsValid():
            shape = checked
        else:
            sf = ShapeFix_Shape(checked)
            sf.Perform()
            if BRepCheck_Analyzer(sf.Shape()).IsValid():
                shape = sf.Shape()

    info["bodies"], info["voids"] = nb, nv
    info["faces"] = count(shape, TopAbs_FACE)
    info["volume"] = volume(shape)
    info["valid"] = BRepCheck_Analyzer(shape).IsValid()
    return shape, info


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
    ap.add_argument("--true-size", action="store_true",
                    help="rebuild at the apparent design size, with radii snapped to round values")
    args = ap.parse_args()

    files = []
    for p in map(Path, args.inputs):
        files += sorted(f for f in p.iterdir() if f.suffix.lower() == ".stl") if p.is_dir() else [p]
    if not files:
        sys.exit("No STL files found.")

    shapes = []
    for f in files:
        t = time.time()
        print(f"{f.name}: converting...", flush=True)
        shape, info = stl_to_solid(f, args.tol, not args.no_fuse, not args.no_curves, args.true_size,
                                   not args.no_blends, not args.no_repair, args.simplify)
        out_dir = Path(args.out) if args.out else f.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / (f.stem + ".step")
        write_step(shape, out)
        shapes.append(shape)
        nb, nv = info["bodies"], info["voids"]
        extra = f" ({nb} bodies{', ' + str(nv) + ' cavities' if nv else ''})" if nb + nv > 1 else ""
        print(f"  {info['triangles']} triangles -> {info['faces']} faces{extra}, "
              f"volume {info['volume']:,.1f} mm^3, "
              f"{'valid solid' if info['valid'] else 'WARNING: check geometry' + (' (the STL itself is not a clean solid: it touches or crosses itself)' if info.get('mesh_defects') else '')}, "
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
