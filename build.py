"""
Build a B-rep shell straight from the mesh and the detected curved patches.

Every flat facet of the mesh becomes one planar face; every detected patch becomes
one exact curved face (swept from its profile, so its outline is the patch's
natural boundary lines). Where a flat face borders a curved patch, that stretch of
its outline is replaced by the exact line or arc, so the faces meet cleanly. All
faces are then sewn together.
"""
import itertools
import math

import numpy as np
from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeVertex)
from OCP.BRepPrimAPI import BRepPrimAPI_MakeRevol, BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeBox
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Splitter
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.GeomAPI import GeomAPI_Interpolate
from OCP.collections import HArray1_gp_Pnt
from OCP.collections import List_TopoDS_Shape
from OCP.GC import GC_MakeArcOfCircle
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopLoc import TopLoc_Location
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.BRepFill import BRepFill_Filling
from OCP.GeomAbs import GeomAbs_C0
from OCP.ShapeFix import ShapeFix_Face, ShapeFix_Shape
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL, TopAbs_REVERSED
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.gp import gp_Ax1, gp_Ax2, gp_Circ, gp_Dir, gp_Pln, gp_Pnt

from features import TWO_PI, _angle_gap, _frame, Revolved, Sphere
from blends import TANGENT_DEG, MAX_DEVIATION, MAX_BULGE, MAX_EDGE_GAP


def _pnt(p):
    return gp_Pnt(float(p[0]), float(p[1]), float(p[2]))


def _line(a, b):
    return BRepBuilderAPI_MakeEdge(_pnt(a), _pnt(b)).Edge()


def _arc(a, mid, b):
    return BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(_pnt(a), _pnt(mid), _pnt(b)).Value()).Edge()


# ---------------------------------------------------------------- where a vertex sits on a patch

class Boundary:
    """Answers "which of the patch's boundary lines is this point on" and builds those lines exactly."""

    def __init__(self, feature, tol):
        self.f, self.m, self.tol = feature, feature.model, tol

    def labels(self, p):
        f, m, tol = self.f, self.m, self.tol
        if f.kind == "wedge":
            return {i for i, n in enumerate(f.planes) if abs((p - m.c) @ n) <= tol}
        if not isinstance(m, Revolved):
            return set()
        rho, z, w = m.local(p)
        rho, z = rho[0], z[0]
        out = set()
        t, scale = self._param(rho, z)
        if abs(t - f.lo) * scale <= tol:
            out.add("lo")
        if abs(t - f.hi) * scale <= tol:
            out.add("hi")
        if f.span < TWO_PI and rho > tol:
            u = m.angle(w)[0]
            if _angle_gap(u, f.u0) * rho <= tol:
                out.add("u0")
            if _angle_gap(u, f.u0 + f.span) * rho <= tol:
                out.add("u1")
        return out

    def _param(self, rho, z):
        m = self.m
        if m.line:
            return z, 1.0
        rc, zc, r = m.circle
        t = math.atan2(z - zc, rho - rc)
        # keep the angle on the same turn as the patch's range
        while t < self.f.lo - math.pi:
            t += TWO_PI
        while t > self.f.hi + math.pi:
            t -= TWO_PI
        return t, r

    def profile_point(self, t, u):
        m = self.m
        if m.line:
            c0, k = m.line
            rho = c0 + k * t
            return m.point(0.0 if abs(rho) < 1e-4 else rho, t, u)
        rc, zc, r = m.circle
        rho = rc + r * math.cos(t)
        return m.point(0.0 if abs(rho) < 1e-4 else rho, zc + r * math.sin(t), u)   # snap onto the axis

    def edge(self, label, pts):
        """Exact edge along boundary line `label` through the run of mesh points `pts`."""
        a, b = pts[0], pts[-1]
        m, f = self.m, self.f
        if f.kind == "wedge":
            mid = pts[len(pts) // 2] if len(pts) > 2 else self._sphere_mid(a, b)
            return _arc(a, mid, b)
        if label in ("u0", "u1"):
            if m.line:
                return _line(a, b)
            mid = pts[len(pts) // 2] if len(pts) > 2 else None
            if mid is None:
                u = f.u0 if label == "u0" else f.u0 + f.span
                ta = self._param(*[x[0] for x in m.local(a)[:2]])[0]
                tb = self._param(*[x[0] for x in m.local(b)[:2]])[0]
                mid = self.profile_point((ta + tb) / 2, u)
            return _arc(a, mid, b)
        # a circle round the axis
        if len(pts) > 2:
            return _arc(a, pts[len(pts) // 2], b)
        ra, za, wa = m.local(a)
        rb, zb, wb = m.local(b)
        ua, ub = m.angle(wa)[0], m.angle(wb)[0]
        du = math.remainder(ub - ua, TWO_PI)
        rho, z = (ra[0] + rb[0]) / 2, (za[0] + zb[0]) / 2
        return _arc(a, m.point(rho, z, ua + du / 2), b)

    def full_circle(self, label):
        m, f = self.m, self.f
        t = f.lo if label == "lo" else f.hi
        p = self.profile_point(t, 0.0)
        rho, z, _ = m.local(p)
        center = m.a + z[0] * m.d
        return BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(_pnt(center), gp_Dir(*m.d)), float(rho[0]))).Edge()

    def _sphere_mid(self, a, b):
        c, r = self.m.c, self.m.r
        v = (a + b) / 2 - c
        return c + v / np.linalg.norm(v) * r


# ---------------------------------------------------------------- mesh boundary loops

def _loops(tris):
    """Directed boundary loops (lists of vertex ids) of a set of consistently wound triangles."""
    directed = {}
    for a, b, c in tris:
        for e in ((a, b), (b, c), (c, a)):
            directed[e] = True
    nxt = {}
    for (a, b) in directed:
        if (b, a) not in directed:
            nxt.setdefault(a, []).append(b)
    loops = []
    while nxt:
        start = next(iter(nxt))
        loop, v = [start], start
        while True:
            outs = nxt.get(v)
            if not outs:
                return None
            w = outs.pop()
            if not outs:
                del nxt[v]
            if w == start:
                break
            loop.append(w)
            v = w
            if len(loop) > len(directed) + 1:
                return None
        loops.append(loop)
    return loops


# ---------------------------------------------------------------- shared edges

class Edges:
    """Builds each boundary run once, so the two faces on either side share one edge.

    A run is a chain of mesh vertices where a face meets one neighbour. It becomes an
    exact line or arc if it follows one of a patch's boundary lines, otherwise a
    smooth spline through the mesh points (which lie on both surfaces).
    """

    def __init__(self, pts):
        self.pts, self.cache = pts, {}

    def get(self, ids, make):
        ids = [int(i) for i in ids]
        if len(ids) > 2 and ids[0] == ids[-1]:
            # a closed loop: canonical start (lowest vertex id) and direction
            ring = ids[:-1]
            i = ring.index(min(ring))
            ring = ring[i:] + ring[:i]
            flip = ring[1] > ring[-1]
            if flip:
                ring = [ring[0]] + ring[1:][::-1]
            canon = ring + [ring[0]]
            key = tuple(canon) + ("closed",)
        else:
            flip = ids[0] > ids[-1]
            canon = ids[::-1] if flip else ids
            key = tuple(canon)
        if key not in self.cache:
            try:
                self.cache[key] = make(canon)
            except Exception:           # e.g. an arc through three nearly aligned points
                self.cache[key] = None
            if self.cache[key] is None and make != self.spline:
                self.cache[key] = self.spline(canon)
        edge = self.cache[key]
        if edge is None:
            return None
        return TopoDS.Edge(edge.Reversed()) if flip else edge

    def spline(self, ids):
        pts = self.pts[ids]
        closed = len(ids) > 3 and ids[0] == ids[-1]
        if closed:
            pts = pts[:-1]
        keep = np.r_[True, np.linalg.norm(np.diff(pts, axis=0), axis=1) > 1e-6]
        pts = pts[keep]
        if len(pts) == 2 and not closed:
            return _line(pts[0], pts[1])
        arr = HArray1_gp_Pnt(1, len(pts))
        for i, p in enumerate(pts):
            arr.SetValue(i + 1, _pnt(p))
        interp = GeomAPI_Interpolate(arr, closed, 1e-7)
        interp.Perform()
        if not interp.IsDone():
            return None
        return BRepBuilderAPI_MakeEdge(interp.Curve()).Edge()


def _run_edge(edges, bounds, tag, ids, patch_side=False):
    """Edge for a run of mesh vertices; tag = (patch, boundary line or None, ...).

    A single mesh edge along a trimmed outline is a chord of the true curve. The flat
    facet next to it keeps the chord; the patch gets an arc bowed onto its own surface
    (sewing closes the tiny gap, which is no bigger than the facet's own error).
    """
    k, label = tag[0], tag[1]
    if label is None and len(ids) == 2:
        if not patch_side or bounds[k].f.kind == "blend":
            return _line(edges.pts[ids[0]], edges.pts[ids[1]])
        m = bounds[k].m
        a, b = edges.pts[ids[0]], edges.pts[ids[1]]
        mid = ((a + b) / 2)[None]
        mid = mid[0] - m.signed(mid)[0] * m.normal(mid)[0]
        return _safe_arc(a, mid, b)
    if label is None:
        return edges.get(ids, edges.spline)
    bound = bounds[k]
    return edges.get(ids, lambda canon: bound.edge(label, edges.pts[canon]))


def _safe_arc(a, mid, b):
    try:
        return _arc(a, mid, b)
    except Exception:
        return _line(a, b)


def _runs(loop, tags):
    """Split a closed loop of vertices into runs of equal tag: [(tag, [vertex ids...]), ...]."""
    n = len(loop)
    start = next((i for i in range(n) if tags[i] != tags[i - 1]), 0)
    out, i = [], 0
    while i < n:
        j = (start + i) % n
        tag = tags[j]
        run = [loop[j]]
        while i < n and tags[(start + i) % n] == tag:
            run.append(loop[(start + i + 1) % n])
            i += 1
            if tag is None:
                break
        out.append((tag, run))
    return out


SHARP_DEG = 30   # a run of mesh vertices turning this sharply at a vertex is split there


def _split_sharp(runs, pts):
    """Split free-form runs (splines through mesh points) at sharp corners: a smooth curve
    forced through a hairpin overshoots wildly. Depends only on the vertex chain, so both
    faces beside a run split it the same way."""
    out = []
    for tag, run in runs:
        if tag is None or tag[1] is not None or len(run) < 3:
            out.append((tag, run))
            continue
        closed = run[0] == run[-1]
        ring = run[:-1] if closed else run
        n = len(ring)

        def sharp(i):
            a, b, c = pts[ring[i - 1]], pts[ring[i]], pts[ring[(i + 1) % n]]
            u, v = b - a, c - b
            nu, nv = np.linalg.norm(u), np.linalg.norm(v)
            return nu > 0 and nv > 0 and u @ v < math.cos(math.radians(SHARP_DEG)) * nu * nv

        corners = [i for i in (range(n) if closed else range(1, n - 1)) if sharp(i)]
        if not corners:
            out.append((tag, run))
            continue
        if closed:
            # start at the corner with the lowest vertex id, so both sides agree
            start = min(corners, key=lambda i: ring[i])
            ring = ring[start:] + ring[:start]
            corners = sorted((i - start) % n for i in corners)
            ring = ring + [ring[0]]
            cuts = corners + [n]
        else:
            ring = list(ring)
            cuts = [0] + corners + [n - 1]
        for i, j in zip(cuts, cuts[1:]):
            out.append((tag, ring[i:j + 1]))
    return out


def _loop_tags(mesh, loop, owner, bounds, edge_tri, patch=None):
    """Tag each edge of a boundary loop.

    None: a plain mesh edge between two flat facets. Otherwise (patch, boundary line):
    a curved patch is involved (the loop's own patch if `patch` is given, else the
    neighbouring one); the boundary line is the patch's line the edge runs along, or
    None if it runs along a trimmed outline. On a patch's own outline a third item
    records the neighbour, so runs split exactly where the neighbouring face's do.
    """
    pts, n, tags = mesh.pts, len(loop), []
    for i in range(n):
        a, b = loop[i], loop[(i + 1) % n]
        other = edge_tri.get((b, a))
        k = owner[mesh.facet_of[other]] if other is not None else -1
        neighbour = ("patch", int(k)) if k >= 0 else (int(mesh.facet_of[other]) if other is not None else None)
        ref = patch if patch is not None else k
        if ref < 0:
            tags.append(None)
            continue
        both = bounds[ref].labels(pts[a]) & bounds[ref].labels(pts[b])
        tags.append((ref, min(both) if both else None, neighbour if patch is not None else None))
    return tags


def _wire(mesh, loop, tags, edges, bounds):
    wire = BRepBuilderAPI_MakeWire()
    pts = mesh.pts
    for tag, run in _split_sharp(_runs(loop, tags), pts):
        if tag is None:
            wire.Add(_line(pts[run[0]], pts[run[1]]))
            continue
        edge = _run_edge(edges, bounds, tag, run)
        if edge is None:
            return None
        wire.Add(edge)
    return wire.Wire() if wire.IsDone() else None


# ---------------------------------------------------------------- faces

def _planar_face(mesh, fid, owner, bounds, edge_tri, edges):
    tris = mesh.tris[mesh.ftris[fid]]
    loops = _loops(tris)
    if loops is None:
        return None
    pts = mesh.pts
    # outer boundary first: the loop enclosing the most area
    normal = mesh.fn[fid]
    area = lambda loop: abs(sum(np.cross(pts[loop[i]], pts[loop[(i + 1) % len(loop)]]) @ normal
                                for i in range(len(loop))))
    loops.sort(key=area, reverse=True)
    wires = []
    for loop in loops:
        n = len(loop)
        tags = _loop_tags(mesh, loop, owner, bounds, edge_tri)
        whole = all(t is not None and t == tags[0] for t in tags)
        if whole and tags[0][1] in ("lo", "hi") and bounds[tags[0][0]].f.kind == "revolve":
            bound = bounds[tags[0][0]]
            edge = bound.full_circle(tags[0][1])
            # the circle runs anticlockwise about the axis; match the loop's direction
            turn = sum(np.cross(pts[loop[i]], pts[loop[(i + 1) % n]]) for i in range(n)) @ bound.m.d
            w = BRepBuilderAPI_MakeWire()
            w.Add(TopoDS.Edge(edge.Reversed()) if turn < 0 else edge)
            wires.append(w.Wire())
            continue
        if whole:
            w = BRepBuilderAPI_MakeWire()
            for tag, run in _split_sharp([(tags[0], loop + [loop[0]])], pts):
                edge = _run_edge(edges, bounds, tag, run)
                if edge is None:
                    return None
                w.Add(edge)
            if not w.IsDone():
                return None
            wires.append(w.Wire())
            continue
        wire = _wire(mesh, loop, tags, edges, bounds)
        if wire is None:
            return None
        wires.append(wire)
    plane = gp_Pln(_pnt(pts[tris[0][0]]), gp_Dir(*normal))
    maker = BRepBuilderAPI_MakeFace(plane, wires[0], True)
    for w in wires[1:]:
        maker.Add(w)
    if not maker.IsDone():
        return None
    fix = ShapeFix_Face(maker.Face())
    fix.Perform()
    return fix.Face()


def _patch_faces(feature, bound, mesh):
    m = feature.model
    if feature.kind == "ball":
        ex = TopExp_Explorer(BRepPrimAPI_MakeSphere(_pnt(m.c), m.r).Shape(), TopAbs_SHELL)
        shell = TopoDS.Shell(ex.Current())
        return [shell if feature.convex else TopoDS.Shell(shell.Reversed())]
    if feature.kind == "wedge":
        return [_wedge_face(feature, bound, mesh)]
    return [_revolved_face(bound, feature.lo, feature.hi, feature.u0, feature.span)]


def _revolved_face(bound, lo, hi, u0, span):
    m = bound.m
    if m.line:
        edge = _line(bound.profile_point(lo, u0), bound.profile_point(hi, u0))
    else:
        edge = _arc(bound.profile_point(lo, u0), bound.profile_point((lo + hi) / 2, u0),
                    bound.profile_point(hi, u0))
    return BRepPrimAPI_MakeRevol(edge, gp_Ax1(_pnt(m.a), gp_Dir(*m.d)), span).Shape()


def _sphere_face(center, r, avoid):
    """Whole sphere face with its seam and poles turned away from the directions `avoid`."""
    avoid = [np.asarray(n, float) for n in avoid]
    mean = np.mean(avoid, axis=0)
    d, e1, e2 = _frame(-mean if np.linalg.norm(mean) > 1e-6 else avoid[0])
    poles = [math.cos(a) * e1 + math.sin(a) * e2 for a in np.radians(np.arange(0, 180, 7.5))]
    pole = max(poles, key=lambda z: min(abs(z @ n) for n in avoid))
    ball = BRepPrimAPI_MakeSphere(gp_Ax2(_pnt(center), gp_Dir(*pole), gp_Dir(*d)), r).Shape()
    return TopoDS.Face(TopExp_Explorer(ball, TopAbs_FACE).Current())


def _generous_surface(feature, bound, mesh):
    """A piece of the patch's surface comfortably bigger than the patch."""
    m = feature.model
    if isinstance(m, Sphere):
        return _sphere_face(m.c, m.r, [mesh.fn[f] for f in feature.facets])
    lo, hi, u0, span = feature.lo, feature.hi, feature.u0, feature.span
    if m.line:
        margin = 0.2 * (hi - lo) + 0.5
        new_lo, new_hi = lo - margin, hi + margin
        c0, k = m.line
        if k != 0:  # don't run a cone past its tip
            apex = -c0 / k
            if new_lo < apex <= lo:
                new_lo = apex
            if hi <= apex < new_hi:
                new_hi = apex
        lo, hi = new_lo, new_hi
    else:
        rc, zc, r = m.circle
        margin = math.radians(10)
        lo, hi = lo - margin, hi + margin
        if hi - lo >= TWO_PI - 1e-6:
            return None
        if (rc + r * np.cos(np.linspace(lo, hi, 60)) < 0).any():
            return None             # would cross the axis
    if span < TWO_PI:
        grow = math.radians(10)
        if span + 2 * grow >= TWO_PI:
            u0, span = 0.0, TWO_PI
        else:
            u0, span = u0 - grow, span + 2 * grow
    return _revolved_face(bound, lo, hi, u0, span)


def _trimmed_face(feature, k, bound, mesh, owner, bounds, edge_tri, edges, tol):
    """A patch cut to an arbitrary outline: split a generous piece of its surface along
    the outline edges and keep the piece the patch's facets lie on."""
    m = feature.model
    tris = np.concatenate([mesh.tris[mesh.ftris[f]] for f in feature.facets])
    loops = _loops(tris)
    if not loops:
        return None
    tools, per_loop = [], []
    for loop in loops:
        tags = _loop_tags(mesh, loop, owner, bounds, edge_tri, patch=k)
        runs = [(tags[0], loop + [loop[0]])] if all(t == tags[0] for t in tags) else _runs(loop, tags)
        per_loop.append([])
        for tag, run in _split_sharp(runs, mesh.pts):
            edge = _run_edge(edges, bounds, tag, run, patch_side=True)
            if edge is None:
                return None
            tools.append(edge)
            per_loop[-1].append(edge)
    base = _generous_surface(feature, bound, mesh)
    if base is None:
        return None
    face = _split_face(feature, mesh, base, tools, loops, tol)
    if face is None:
        face = _outline_face(feature, mesh, base, loops, per_loop)
    return face if face is not None else _polygon_face(feature, mesh, base, loops)


def _outline_face(feature, mesh, base, loops, per_loop):
    """Fallback: the patch's own outline edges (shared with its neighbours) laid on its
    surface as wires, the face fitted to them by ShapeFix."""
    surface = BRep_Tool.Surface_s(TopoDS.Face(TopExp_Explorer(base, TopAbs_FACE).Current()))
    target = float(mesh.farea[feature.facets].sum())
    normal = mesh.fn[feature.facets].T @ mesh.farea[feature.facets]
    pts = mesh.pts
    area = lambda loop: abs(sum(np.cross(pts[loop[i]], pts[loop[(i + 1) % len(loop)]]) @ normal
                                for i in range(len(loop))))
    order = sorted(range(len(loops)), key=lambda i: -area(loops[i]))
    wires = []
    for i in order:
        w = BRepBuilderAPI_MakeWire()
        for e in per_loop[i]:
            w.Add(e)
        if not w.IsDone():
            return None
        wires.append(w.Wire())
    for flip in (False, True):
        ws = [TopoDS.Wire(w.Reversed()) for w in wires] if flip else wires
        maker = BRepBuilderAPI_MakeFace(surface, ws[0], False)
        for w in ws[1:]:
            maker.Add(w)
        if not maker.IsDone():
            continue
        fix = ShapeFix_Face(maker.Face())
        fix.Perform()
        face = fix.Face()
        if not BRepCheck_Analyzer(face).IsValid():
            continue
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props)
        if abs(props.Mass() - target) <= 0.05 * target + 1e-3:
            return face
    return None


def _polygon_face(feature, mesh, base, loops):
    """Fallback (the stlToSolid project's method): lay the patch's own mesh outline,
    as straight edges, on its surface and let ShapeFix fit the face to it."""
    surface = BRep_Tool.Surface_s(TopoDS.Face(TopExp_Explorer(base, TopAbs_FACE).Current()))
    target = float(mesh.farea[feature.facets].sum())
    loops = sorted(loops, key=len, reverse=True)
    for flip in (False, True):
        wires = []
        for loop in loops:
            ring = loop[::-1] if flip else loop
            w = BRepBuilderAPI_MakeWire()
            for a, b in zip(ring, ring[1:] + ring[:1]):
                w.Add(_line(mesh.pts[a], mesh.pts[b]))
            if not w.IsDone():
                return None
            wires.append(w.Wire())
        maker = BRepBuilderAPI_MakeFace(surface, wires[0], False)
        for w in wires[1:]:
            maker.Add(w)
        if not maker.IsDone():
            continue
        fix = ShapeFix_Face(maker.Face())
        fix.Perform()
        face = fix.Face()
        if not BRepCheck_Analyzer(face).IsValid():
            continue
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props)
        if abs(props.Mass() - target) <= 0.05 * target + 1e-3:
            return face
    return None


def _split_face(feature, mesh, base, tools, loops, tol):
    m = feature.model
    splitter = BRepAlgoAPI_Splitter()
    splitter.SetArguments(_shapes([base]))
    splitter.SetTools(_shapes(tools))
    splitter.SetFuzzyValue(max(10 * tol, 1e-4))
    splitter.Build()
    if not splitter.IsDone():
        return None

    # Which pieces make up the patch: drop points from the patch's facets onto the true
    # surface and keep every piece they land on (a patch can come out as several pieces).
    facets = feature.facets
    # sample away from the outline: a needle-thin facet on the edge touches both pieces
    rim = set(np.concatenate(loops).tolist())
    inner = np.array([f for f in facets if not rim & set(mesh.fverts[f].tolist())])
    source = inner if len(inner) else facets
    pick = source[np.linspace(0, len(source) - 1, min(len(source), 60)).astype(int)]
    c = mesh.fcent[pick]
    probes = c - m.signed(c)[:, None] * m.normal(c)
    vertices = [BRepBuilderAPI_MakeVertex(_pnt(p)).Vertex() for p in probes]
    hit = []
    ex = TopExp_Explorer(splitter.Shape(), TopAbs_FACE)
    while ex.More():
        face = TopoDS.Face(ex.Current())
        ex.Next()
        if any(BRepExtrema_DistShapeShape(v, face).Value() < 1e-3 for v in vertices):
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, props)
            hit.append((face, props.Mass()))
    target = float(mesh.farea[facets].sum())
    if not hit or len(hit) > 10:
        return None
    # Normally every piece hit belongs to the patch; if they don't add up (a probe sat
    # right on the outline, touching the next piece too), take the best-matching set.
    total = lambda combo: sum(a for _, a in combo)
    if abs(total(hit) - target) <= 0.1 * target:
        kept = [f for f, _ in hit]
    else:
        combos = [c for n in range(1, len(hit) + 1) for c in itertools.combinations(hit, n)
                  if abs(total(c) - target) <= 0.1 * target]
        best = min(combos, key=lambda c: abs(total(c) - target), default=())
        kept = [f for f, _ in best]
    if not kept:
        return None
    if len(kept) == 1:
        return kept[0]
    comp = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(comp)
    for face in kept:
        builder.Add(comp, face)
    return comp


def _repaired(face):
    """The face if valid, else the first standard repair that makes it valid (edge
    tolerances widened to the real gap, at most twice the allowed deviation)."""
    if BRepCheck_Analyzer(face).IsValid():
        return face
    fix = ShapeFix_Face(face)
    fix.Perform()
    if BRepCheck_Analyzer(fix.Face()).IsValid():
        return fix.Face()
    fix = ShapeFix_Shape(face)
    fix.SetMaxTolerance(2 * MAX_DEVIATION)
    fix.Perform()
    ex = TopExp_Explorer(fix.Shape(), TopAbs_FACE)
    if not ex.More():
        return None
    out = TopoDS.Face(ex.Current())
    ex.Next()
    return out if not ex.More() and BRepCheck_Analyzer(out).IsValid() else None


# degree, points per boundary curve, iterations, anisotropy, 2d/3d/angular/curvature
# tolerances, max degree, max segments
USE_GUIDES = False   # guide points past the outline pull against the real points where the neighbour is curved
BULGE_FACTOR = 1.0   # how far past a circular arc's sag a blend may bow over a facet
FILL_SETTINGS = (3, 15, 2, False, 1e-5, 1e-4, 1e-2, 0.1, 8, 9)
FINE_FILL_SETTINGS = (3, 30, 3, False, 1e-5, 1e-5, 1e-2, 0.1, 8, 20)
FINE_BLENDS = 150       # a part with more blends than this is mostly freeform: skip the finer (slower) refits


def _blend_face(feature, k, mesh, owner, bounds, edge_tri, edges):
    """A smooth freeform face for a blend: an N-sided patch spanning the blend's outline
    (the same edges its neighbours use) through the mesh corners inside it, plus a few
    guide points just past the outline so it rolls tangentially into the faces beside
    it. Sets the feature's volume change and fit on the way."""
    facets = feature.facets
    inside = set(facets.tolist())
    tids = np.concatenate([mesh.ftris[f] for f in facets])
    tris = mesh.tris[tids]
    loops = _loops(tris)
    if not loops or len(loops) != 1:
        return None                 # one outline only (a ring-shaped area is left as it was)
    loop = loops[0]
    pts = mesh.pts
    V = np.unique(tris)
    rim = set(loop)
    # typical facet width, to place the guide points
    ab = np.linalg.norm(pts[tris[:, 1]] - pts[tris[:, 0]], axis=1)
    width = float(np.median(2 * mesh.tarea[tids] / np.maximum(ab, 1e-9)))
    guides = []
    for a, b in zip(loop, loop[1:] + loop[:1]):
        t_out, t_in = edge_tri.get((b, a)), edge_tri.get((a, b))
        if t_out is None or t_in is None:
            continue
        f_in, f_out = mesh.facet_of[t_in], mesh.facet_of[t_out]
        if mesh.fn[f_in] @ mesh.fn[f_out] < math.cos(math.radians(TANGENT_DEG)):
            continue                # a crease: nothing to roll into
        e = pts[b] - pts[a]
        e = e / np.linalg.norm(e)
        mid = (pts[a] + pts[b]) / 2
        away = np.cross(mesh.fn[f_out], e)
        if away @ (pts[mesh.tris[t_in]].mean(axis=0) - mid) > 0:
            away = -away
        guides.append(mid + 0.35 * max(width, 0.05) * away)

    tags = _loop_tags(mesh, loop, owner, bounds, edge_tri, patch=k)
    runs = [(tags[0], loop + [loop[0]])] if all(t == tags[0] for t in tags) else _runs(loop, tags)
    runs = _split_sharp(runs, pts)
    # The same outline gives the same face: reuse it from an earlier attempt at the part
    cache = mesh.__dict__.setdefault("blend_cache", {})
    key = (facets.tobytes(), tuple(tuple(int(v) for v in run) for _, run in runs))
    if key in cache:
        if cache[key] is None:
            return None
        face, feature.change, feature.tolerance, feature.worst, feature.detail = cache[key]
        return face
    cache[key] = None
    outline = []
    for tag, run in runs:
        edge = _run_edge(edges, bounds, tag, run, patch_side=True)
        if edge is None:
            return None
        outline.append(edge)
    inner = [pts[v] for v in V if v not in rim] + (guides if USE_GUIDES else [])

    def filled(settings):
        fill = BRepFill_Filling(*settings)
        for edge in outline:
            fill.Add(edge, GeomAbs_C0, True)
        for q in inner:
            fill.Add(_pnt(q))
        try:
            fill.Build()
            return _repaired(fill.Face()) if fill.IsDone() else None
        except Exception:
            return None

    # a quick patch first; if it follows its outline loosely, a finer (slower) one
    face = filled(FILL_SETTINGS)
    if (face is None or _edge_gap(face) > MAX_EDGE_GAP / 4) and mesh.__dict__.get("blend_count", 0) <= FINE_BLENDS:
        finer = filled(FINE_FILL_SETTINGS)
        if finer is not None and (face is None or _edge_gap(finer) < _edge_gap(face)):
            face = finer
    if face is None or _edge_gap(face) > MAX_EDGE_GAP:
        return None     # a loose outline would force a loose sewing tolerance on the whole part
    surface = BRep_Tool.Surface_s(face)
    proj = GeomAPI_ProjectPointOnSurf()
    proj.Init(surface, *surface.Bounds())

    def along(p, n):
        """Signed distance from p to the surface, measured along direction n."""
        proj.Perform(_pnt(p))
        if not proj.IsDone() or not proj.NbPoints():
            return None
        q = proj.NearestPoint()
        return float((np.array([q.X(), q.Y(), q.Z()]) - p) @ n)

    # how far the surface strays from the mesh corners, and the volume it adds
    worst, change, spread = 0.0, 0.0, 0.0
    for t in tids:
        n_t = mesh.tn[t]
        corners = pts[mesh.tris[t]]
        for p in corners:
            d = along(p, n_t)
            if d is None:
                return None
            worst = max(worst, abs(d))
        s_mid = [along(p, n_t) for p in (corners + corners[[1, 2, 0]]) / 2]
        # between the corners a smooth surface bows away from a flat facet only a little
        # (about a tenth of the facet's size even for a fillet cut into two strips);
        # more than that is the fit overshooting
        centre = corners.mean(axis=0)
        s_in = [along(p, n_t) for p in [centre] + list((corners + centre) / 2)]
        if any(d is None for d in s_mid + s_in):
            return None
        size = float(np.linalg.norm(corners - corners[[1, 2, 0]], axis=1).max())
        if max(abs(d) for d in s_mid + s_in) > max(MAX_DEVIATION, MAX_BULGE * size):
            return None
        change += mesh.tarea[t] * float(np.mean(s_mid))
        spread += mesh.tarea[t] * float(np.mean(np.abs(s_mid)))
    if worst > MAX_DEVIATION:
        return None
    if not _hugs_mesh(face, mesh, tids):
        return None
    gp = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, gp)
    target = float(mesh.farea[facets].sum())
    if abs(gp.Mass() - target) > 0.15 * target + 1e-3:
        return None
    # face the same way as the mesh
    proj.Perform(_pnt(mesh.fcent[facets[0]]))
    u, v = proj.LowerDistanceParameters()
    props = GeomLProp_SLProps(surface, u, v, 1, 1e-9)
    if props.IsNormalDefined():
        n = props.Normal()
        n = np.array([n.X(), n.Y(), n.Z()]) * (-1 if face.Orientation() == TopAbs_REVERSED else 1)
        if n @ mesh.fn[facets[0]] < 0:
            face = TopoDS.Face(face.Reversed())
    feature.change = change
    feature.tolerance = 0.5 * spread + 1e-4 * target
    # its outline is shared with its neighbours, so the gap sewing must close is just
    # how far the face's own edges stray (not how far it lies from the mesh corners)
    feature.worst = _edge_gap(face)
    feature.detail = feature.detail.split(", within")[0] + f", within {worst:.3f} mm of the mesh"
    cache[key] = (face, feature.change, feature.tolerance, feature.worst, feature.detail)
    return face


def _point_triangle_distance(P, T):
    """Distance from each point in P (n, 3) to the nearest of triangles T (m, 3, 3), and
    which triangle that is."""
    best = np.full(len(P), np.inf)
    which = np.zeros(len(P), int)
    for start in range(0, len(T), 256):
        tri = T[start:start + 256]
        a, b, c = tri[:, 0][None], tri[:, 1][None], tri[:, 2][None]
        p = P[:, None]
        ab, ac, ap = b - a, c - a, p - a
        d1, d2 = (ab * ap).sum(-1), (ac * ap).sum(-1)
        bp = p - b
        d3, d4 = (ab * bp).sum(-1), (ac * bp).sum(-1)
        cp = p - c
        d5, d6 = (ab * cp).sum(-1), (ac * cp).sum(-1)
        va = d3 * d6 - d5 * d4
        vb = d5 * d2 - d1 * d6
        vc = d1 * d4 - d3 * d2
        denom = np.where(np.abs(va + vb + vc) > 1e-30, va + vb + vc, 1e-30)
        v, w = vb / denom, vc / denom
        inside = (v >= 0) & (w >= 0) & (v + w <= 1)
        # inside: straight down onto the triangle; outside: the nearest point of its edges
        dist = np.where(inside, np.linalg.norm(p - (a + ab * v[..., None] + ac * w[..., None]), axis=-1), np.inf)
        for e0, e1 in ((a, b), (b, c), (c, a)):
            e = e1 - e0
            t = np.clip(((p - e0) * e).sum(-1) / np.maximum((e * e).sum(-1), 1e-30), 0, 1)
            edge = np.linalg.norm(p - (e0 + e * t[..., None]), axis=-1)
            dist = np.where(inside, dist, np.minimum(dist, edge))
        k = dist.argmin(axis=1)
        dmin = dist[np.arange(len(P)), k]
        better = dmin < best
        best[better], which[better] = dmin[better], k[better] + start
    return best, which


def _hugs_mesh(face, mesh, tids):
    """Does the finished face stay as close to the mesh as a smooth surface through its
    corners can? A smooth surface bows away from a flat facet by only about L*theta/4
    (L its size, theta how far the surface turns across it, read from the mesh's own
    corner normals); a fitted patch that waves between the corners bows much further."""
    BRepMesh_IncrementalMesh(face, 0.005, False, 0.2, True)
    loc = TopLoc_Location()
    tri = BRep_Tool.Triangulation_s(face, loc)
    if tri is None:
        return False
    nodes = np.array([[tri.Node(i).X(), tri.Node(i).Y(), tri.Node(i).Z()] for i in range(1, tri.NbNodes() + 1)])
    ft = np.array([[tri.Triangle(i).Value(j) for j in (1, 2, 3)] for i in range(1, tri.NbTriangles() + 1)]) - 1
    samples = np.vstack([nodes, nodes[ft].mean(axis=1)])
    samples = samples[::max(1, len(samples) // 4000)]
    T = mesh.pts[mesh.tris[tids]]
    d, k = _point_triangle_distance(samples, T)
    # corner normals from the blend's own triangles only (averaging in a face across a
    # sharp edge would make the surface look far more curved than it is)
    corners = mesh.tris[tids]
    vn = np.zeros_like(mesh.pts)
    for j in range(3):
        np.add.at(vn, corners[:, j], mesh.tn[tids] * mesh.tarea[tids][:, None])
    vn /= np.maximum(np.linalg.norm(vn, axis=1), 1e-12)[:, None]
    theta = np.arccos(np.clip(np.einsum("tkj,tj->tk", vn[corners], mesh.tn[tids]), -1, 1)).max(axis=1)
    size = np.linalg.norm(T - T[:, [1, 2, 0]], axis=2).max(axis=1)
    allowed = BULGE_FACTOR * size * theta / 4 + 0.005
    if (d > allowed[k]).any():
        return False
    # and every mesh corner lies on the face itself (not just on its untrimmed surface)
    dc, _ = _point_triangle_distance(mesh.pts[np.unique(corners)], nodes[ft])
    return bool((dc <= MAX_DEVIATION).all())


def _edge_gap(face):
    gap, ex = 0.0, TopExp_Explorer(face, TopAbs_EDGE)
    while ex.More():
        gap = max(gap, BRep_Tool.Tolerance_s(TopoDS.Edge(ex.Current())))
        ex.Next()
    return gap


def _shapes(items):
    out = List_TopoDS_Shape()
    for x in items:
        out.Append(x)
    return out


def _wedge_face(feature, bound, mesh):
    """A sphere face trimmed by the planes through its centre that bound the corner."""
    m = feature.model
    mean = np.mean([np.asarray(n) for n in feature.planes], axis=0)
    d, e1, e2 = _frame(-mean if np.linalg.norm(mean) > 1e-6 else feature.planes[0])
    # Keep the sphere's seam (on the far side, along d) and its poles well away from
    # the corner patch and from every trimming plane: a plane through a pole fails.
    poles = [math.cos(a) * e1 + math.sin(a) * e2 for a in np.radians(np.arange(0, 180, 7.5))]
    pole = max(poles, key=lambda z: min(abs(z @ n) for n in feature.planes))
    ball = BRepPrimAPI_MakeSphere(gp_Ax2(_pnt(m.c), gp_Dir(*pole), gp_Dir(*d)), m.r).Shape()
    ex = TopExp_Explorer(ball, TopAbs_FACE)
    shape = TopoDS.Face(ex.Current())
    size = 2 * m.r + 1
    for n in feature.planes:
        nd, n1, n2 = _frame(n)
        corner = m.c - size * n1 - size * n2
        box = BRepPrimAPI_MakeBox(gp_Ax2(_pnt(corner), gp_Dir(*nd), gp_Dir(*n1)), 2 * size, 2 * size, size).Shape()
        common = BRepAlgoAPI_Common(shape, box)
        if not common.IsDone():
            return None
        shape = common.Shape()
    return shape


def settle_blends(mesh, features, tol, split, skipped):
    """Build every blend's face on its own first (each depends only on its own outline),
    cutting in two any that won't fit and giving back the pieces of any that can't be
    cut further, so the whole part is then built and sewn just once or twice."""
    features = list(features)
    mesh.blend_count = sum(f.kind == "blend" for f in features)
    for _ in range(4):
        owner = np.full(len(mesh.fn), -1)
        for k, f in enumerate(features):
            owner[f.facets] = k
        bounds = [Boundary(f, 10 * tol) for f in features]
        edge_tri = {}
        for t, (a, b, c) in enumerate(mesh.tris):
            edge_tri[(a, b)] = t
            edge_tri[(b, c)] = t
            edge_tri[(c, a)] = t
        edges = Edges(mesh.pts)
        out, changed = [], False
        for k, f in enumerate(features):
            if f.kind != "blend" or _blend_face(f, k, mesh, owner, bounds, edge_tri, edges) is not None:
                out.append(f)
                continue
            changed = True
            halves = split(mesh, f)
            if halves:
                out += halves
            else:
                skipped.append(f)
                out += list(f.parts)
        features = out
        if not changed:
            break
    return features


def build_faces(mesh, features, tol):
    """All faces of the rebuilt part, as a compound, plus the features that were used."""
    owner = np.full(len(mesh.fn), -1)
    for k, f in enumerate(features):
        owner[f.facets] = k
    bounds = [Boundary(f, 10 * tol) for f in features]
    edge_tri = {}
    for t, (a, b, c) in enumerate(mesh.tris):
        edge_tri[(a, b)] = t
        edge_tri[(b, c)] = t
        edge_tri[(c, a)] = t

    comp = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(comp)
    failed, shells = [], []
    edges = Edges(mesh.pts)
    for k, (f, bound) in enumerate(zip(features, bounds)):
        if f.kind == "trimmed":
            faces = [_trimmed_face(f, k, bound, mesh, owner, bounds, edge_tri, edges, tol)]
        elif f.kind == "blend":
            faces = [_blend_face(f, k, mesh, owner, bounds, edge_tri, edges)]
        else:
            faces = _patch_faces(f, bound, mesh)
        if any(x is None for x in faces):
            failed.append(k)
            continue
        for x in faces:
            if f.kind == "ball":
                shells.append(x)
            else:
                builder.Add(comp, x)
    for fid in range(len(mesh.fn)):
        if owner[fid] >= 0 and owner[fid] not in failed:
            continue
        face = None if owner[fid] >= 0 else _planar_face(mesh, fid, owner, bounds, edge_tri, edges)
        if face is None:
            # fall back to the facet's own triangles
            for t in mesh.ftris[fid]:
                a, b, c = mesh.pts[mesh.tris[t]]
                w = BRepBuilderAPI_MakeWire(_line(a, b), _line(b, c), _line(c, a)).Wire()
                builder.Add(comp, BRepBuilderAPI_MakeFace(w, True).Face())
        else:
            builder.Add(comp, face)
    return comp, shells, failed


def sew(shape, tol, shells=()):
    """Sew the faces into shells; closed shells made elsewhere (whole spheres) are added as
    they are. Returns (shape, midpoints of the edges left unmatched)."""
    s = BRepBuilderAPI_Sewing(tol)
    s.Add(shape)
    s.Perform()
    free = []
    for i in range(1, s.NbFreeEdges() + 1):
        c = BRepAdaptor_Curve(s.FreeEdge(i))
        p = c.Value((c.FirstParameter() + c.LastParameter()) / 2)
        free.append((p.X(), p.Y(), p.Z()))
    if not shells:
        return s.SewedShape(), free
    comp = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(comp)
    builder.Add(comp, s.SewedShape())
    for shell in shells:
        builder.Add(comp, shell)
    return comp, free


def features_near(mesh, features, points, reach=0.5):
    """Indices of the features whose facets come within `reach` mm of any of the points."""
    if not points:
        return []
    points = np.asarray(points)
    out = []
    for k, f in enumerate(features):
        P = mesh.pts[np.unique(np.concatenate([mesh.fverts[x] for x in f.facets]))]
        lo, hi = P.min(axis=0) - reach, P.max(axis=0) + reach
        near = points[np.all((points >= lo) & (points <= hi), axis=1)]
        if len(near) and (np.linalg.norm(near[:, None] - P[None], axis=2).min() <= reach):
            out.append(k)
    return out
