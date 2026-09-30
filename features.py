"""
Find curved surfaces in a triangle mesh so they can be rebuilt as true curves.

The mesh is split into patches, each lying on one simple surface:

  cylinder  round holes, pins, straight rounded edges (fillets)
  cone      countersinks, chamfers around round edges, cone tips
  sphere    domes, dimples, balls, rounded corners where three edges meet
  torus     rounded edges that follow a curve (the corners of a rounded box lid,
            a fillet around the base of a pin, a rounded rim on a hole)

How a patch is found: take two or three neighbouring facets, work out the one
surface they could lie on, then grow outwards while the facets' corners stay on
that surface. Straight chamfers and flat faces need nothing: they are already
exact planes.

A patch is only accepted if its outline follows the surface's natural boundary
lines (for a cylinder: the two end circles and the two straight sides). This
rejects partial matches, holes crossed by other features, oblique cuts and so on;
those areas simply stay faceted.

Each accepted patch records its exact surface and its boundary lines; build.py
turns it into one true curved face and stitches it to the flat faces around it.
"""
import math
import struct
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np

TOL = 0.001          # mm: CAD exports put mesh corners on the true surface to within this
FLOAT_TOL = 2e-7     # plus this fraction of the largest coordinate (STL stores 32-bit floats)
NORMAL_DEG = 30      # a facet may tilt this far from the true surface normal (coarse meshes)
SMOOTH_DEG = 40      # facets meeting at a gentler bend than this can be one curved surface
COPLANAR_DEG = 0.05  # triangles closer to flat than this are merged into one facet
MAX_SAG = 0.25       # mm: the middle of a facet may sit this far inside the curve (coarse meshes)
BEND_SLACK_DEG = 8   # neighbouring facets must bend by what the surface predicts, within this
MIN_SPAN_DEG = 15    # a partial patch must curve through at least this much
CREASE_DEG = 20      # a one-ring cone must meet a neighbour at a crease at least this sharp
MIN_CORNERS = {"cylinder": 8, "cone": 10, "sphere": 8, "torus": 12}
# A patch whose outline doesn't follow its own boundary lines (e.g. a hole running
# out through a sloped face, or a rounded knob cut by another feature) is accepted as
# a "trimmed" patch. Without the outline check to lean on, it needs more evidence:
WHOLE_ON_SURFACE = 0.6     # a region fitted whole: this share of its corners on the surface
WHOLE_STRAY = 0.02          # mm: how far the rest (outline corners on a seam) may stray
TRIMMED_MIN_CORNERS = 12
TRIMMED_MIN_FACETS = 6
TRIMMED_MIN_TURN_DEG = 15   # the facets must face in directions at least this far apart

# Loose pass (run last, on leftovers only; the approach of the stlToSolid project):
# accept any patch the surface explains geometrically, without the outline and crease
# rules. Accuracy is still guaranteed by the corner, interior and normal checks; the
# strict passes run first so clean patches keep their exact natural boundaries.
LOOSE_NORMAL_DEG = 8        # facet normal vs surface normal (plus the facet's own spread)
LOOSE_MAX_GAP = 0.08        # mm: facet interiors (centres, edge midpoints) vs surface
LOOSE_MIN_TURN_DEG = {"cylinder": 6, "cone": 6, "sphere": 20, "torus": 20}
LOOSE_MIN_FACETS = {"cylinder": 3, "cone": 6, "sphere": 8, "torus": 8}
_loose = False
_anchored = False   # the surface is pinned by the faces around it (fillets.py): trust it
_axis_vouched = False   # on an axis another patch already has, in a rounded-off mesh: a cut-to-
                        # shape patch may hand over smoothly (to its own fillets) all round
TWO_PI = 2 * math.pi


# ---------------------------------------------------------------- mesh

def load_stl(path):
    """Return welded vertex array and triangle index array for a binary or ASCII STL."""
    data = Path(path).read_bytes()
    n = struct.unpack("<I", data[80:84])[0] if len(data) >= 84 else -1
    if len(data) == 84 + 50 * n:
        rec = np.frombuffer(data, offset=84, count=n, dtype=np.dtype(
            [("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
        verts = rec["v"].reshape(-1, 3).astype(float)
    else:
        lines = data.decode(errors="ignore").splitlines()
        verts = np.array([l.split()[1:4] for l in lines if l.strip().startswith("vertex")], float)
    key = np.round(verts * 1e4).astype(np.int64)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    pts = np.zeros((len(uniq), 3))
    pts[inv.ravel()] = verts
    tris = inv.reshape(-1, 3)
    tris = tris[(tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 2] != tris[:, 0])]
    return pts, tris


def _labels(n, a, b):
    """Connected-component label for each of n items linked by pairs (a[i], b[i])."""
    parent = np.arange(n)

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in zip(a, b):
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
    roots = np.array([root(i) for i in range(n)])
    return np.unique(roots, return_inverse=True)[1]


def _groups(labels):
    order = np.argsort(labels, kind="stable")
    return np.split(order, np.nonzero(np.diff(labels[order]))[0] + 1)


def _edges(tris):
    return np.sort(tris[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)


MERGE_MAX_DEG = 1     # rounding noise never explains a bigger tilt than this
AXIS_NOISE = 10       # on a known axis, corners may sit this many rounding steps off
MIN_VOUCHED_CONE_DEG = 5    # ... but a cone found that way must taper at least this much


def grid_noise(pts):
    """How far corners may sit off the true surface because the file rounded their
    coordinates to a grid (many exports write 0.001 mm steps): half a step along each
    axis. 0 when the coordinates aren't rounded."""
    for step in (0.01, 0.001, 0.0001):
        err = np.abs(pts - np.round(pts / step) * step)
        if (err <= 2.5e-7 * np.abs(pts) + 1e-7).mean() >= 0.75:
            return step * math.sqrt(3) / 2
    return 0.0


def _merge_flat(pts, tris, tn, tarea, labels, t1, t2, noise, tol):
    """Join neighbouring facets that differ only by rounding noise. Rounded corners tilt
    a long thin triangle by up to noise / its width, which splits one flat face into many
    facets. Two triangles are joined when their tilt is within what the noise explains,
    and only while all the corners of the joined facet still lie on one plane (so a
    finely cut curve can't creep into a plane one strip at a time)."""
    if noise <= 0 or not len(t1):
        return labels
    corners = pts[tris]
    longest = np.linalg.norm(corners - corners[:, [1, 2, 0]], axis=2).max(axis=1)
    width = 2 * tarea / np.maximum(longest, 1e-12)
    dot = np.clip(np.einsum("ij,ij->i", tn[t1], tn[t2]), -1, 1)
    allowed = np.minimum(noise / np.maximum(width[t1], 1e-12) + noise / np.maximum(width[t2], 1e-12),
                         math.radians(MERGE_MAX_DEG))    # (specks narrower than the noise: no)
    ok = (np.arccos(dot) <= allowed) & (labels[t1] != labels[t2])
    if not ok.any():
        return labels
    verts = {}
    for t, f in enumerate(labels):
        verts.setdefault(int(f), set()).update(tris[t].tolist())
    parent = list(range(labels.max() + 1))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in np.nonzero(ok)[0][np.argsort(-dot[ok])]:
        x, y = root(int(labels[t1[i]])), root(int(labels[t2[i]]))
        if x == y:
            continue
        both = verts[x] | verts[y]
        P = pts[list(both)]
        mid = P.mean(axis=0)
        n = np.linalg.svd(P - mid, full_matrices=False)[2][2]
        if np.abs((P - mid) @ n).max() <= tol + noise:
            parent[x] = y
            verts[y] = both
            del verts[x]
    return np.unique([root(int(f)) for f in labels], return_inverse=True)[1]


class Mesh:
    """Triangles grouped into flat facets, with which facets meet at a gentle bend."""

    def __init__(self, pts, tris):
        a, b, c = (pts[tris[:, k]] for k in range(3))
        cross = np.cross(b - a, c - a)
        size = np.linalg.norm(cross, axis=1)
        keep = size > 1e-12
        self.pts, self.tris = pts, tris[keep]
        self.tn = cross[keep] / size[keep, None]
        self.tarea = size[keep] / 2

        edges = _edges(self.tris)
        owner = np.repeat(np.arange(len(self.tris)), 3)
        order = np.lexsort((edges[:, 1], edges[:, 0]))
        edges, owner = edges[order], owner[order]
        shared = np.nonzero(np.all(edges[1:] == edges[:-1], axis=1))[0]
        t1, t2 = owner[shared], owner[shared + 1]
        dot = np.einsum("ij,ij->i", self.tn[t1], self.tn[t2])
        flat = dot > math.cos(math.radians(COPLANAR_DEG))

        self.facet_of = _labels(len(self.tris), t1[flat], t2[flat])
        self.noise = grid_noise(pts)
        self.facet_of = _merge_flat(self.pts, self.tris, self.tn, self.tarea, self.facet_of,
                                    t1[~flat], t2[~flat], self.noise, TOL + FLOAT_TOL * float(np.abs(pts).max()))
        nf = self.facet_of.max() + 1
        self.ftris = _groups(self.facet_of)
        self.farea = np.bincount(self.facet_of, self.tarea, nf)
        fn = np.zeros((nf, 3))
        np.add.at(fn, self.facet_of, self.tn * self.tarea[:, None])
        self.fn = fn / np.linalg.norm(fn, axis=1)[:, None]
        fc = np.zeros((nf, 3))
        np.add.at(fc, self.facet_of, self.pts[self.tris].mean(axis=1) * self.tarea[:, None])
        self.fcent = fc / self.farea[:, None]
        self.fverts = [np.unique(self.tris[t]) for t in self.ftris]

        f1, f2 = self.facet_of[t1[~flat]], self.facet_of[t2[~flat]]
        ends = edges[shared[~flat]]
        pairs = np.sort(np.c_[f1, f2], axis=1)
        bend = np.einsum("ij,ij->i", self.fn[pairs[:, 0]], self.fn[pairs[:, 1]])
        smooth = (pairs[:, 0] != pairs[:, 1]) & (bend > math.cos(math.radians(SMOOTH_DEG)))
        self.edge_facets = {}   # (vertex, vertex), sorted -> facets using that edge
        for (a, b), f in zip(edges, self.facet_of[owner]):
            self.edge_facets.setdefault((a, b), []).append(f)
        self.nbrs = [set() for _ in range(nf)]
        self.shared = {}   # (facet, facet) -> list of shared edges (vertex id pairs)
        for (x, y), e in zip(pairs[smooth], ends[smooth]):
            self.nbrs[x].add(y)
            self.nbrs[y].add(x)
            self.shared.setdefault((x, y), []).append(e)

    def smooth_regions(self):
        nf = len(self.fn)
        a = [x for x in range(nf) for y in self.nbrs[x] if x < y]
        b = [y for x in range(nf) for y in self.nbrs[x] if x < y]
        return [g for g in _groups(_labels(nf, a, b)) if len(g) >= 3]

    def shared_direction(self, f, g):
        """Direction of the (longest) edge shared by two facets."""
        e = self.shared.get((min(f, g), max(f, g)))
        if not e:
            return None
        vec = [self.pts[j] - self.pts[i] for i, j in e]
        v = max(vec, key=np.linalg.norm)
        return v / np.linalg.norm(v)


# ---------------------------------------------------------------- surfaces

def _frame(d):
    d = np.asarray(d, float) / np.linalg.norm(d)
    e1 = np.cross(d, [1.0, 0, 0] if abs(d[0]) < 0.9 else [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    return d, e1, np.cross(d, e1)


class Revolved:
    """Surface of revolution: cylinder, cone, sphere or torus.

    Around the axis (point a, direction d), a point has radius rho and height z.
    The profile is a line rho = c0 + k*z (cylinder when k == 0, else cone) or a
    circle of radius r centred at (rc, zc) (sphere when rc == 0, else torus).
    """

    def __init__(self, a, d, line=None, circle=None):
        self.d, self.e1, self.e2 = _frame(d)
        a = np.asarray(a, float)
        self.a = a - (a @ self.d) * self.d
        self.line, self.circle = line, circle
        if line:
            self.scale = max(abs(line[0]), 1.0)
        else:
            self.scale = circle[2]

    @property
    def kind(self):
        if self.line:
            return "cylinder" if self.line[1] == 0 else "cone"
        return "sphere" if self.circle[0] == 0 else "torus"

    def local(self, p):
        v = np.atleast_2d(p) - self.a
        z = v @ self.d
        w = v - np.outer(z, self.d)
        return np.linalg.norm(w, axis=1), z, w

    def angle(self, w):
        return np.arctan2(w @ self.e2, w @ self.e1)

    def signed(self, p):
        """Distance to the surface; positive = away from the axis (line) / tube centre (circle)."""
        rho, z, _ = self.local(p)
        if self.line:
            c0, k = self.line
            return (rho - c0 - k * z) / math.hypot(1, k)
        rc, zc, r = self.circle
        return np.hypot(rho - rc, z - zc) - r

    def normal(self, p):
        rho, z, w = self.local(p)
        radial = w / np.maximum(rho, 1e-12)[:, None]
        radial[rho < 1e-9] = 0
        if self.line:
            c0, k = self.line
            nr = np.full(len(rho), 1 / math.hypot(1, k))
            nz = np.full(len(rho), -k / math.hypot(1, k))
        else:
            rc, zc, r = self.circle
            nr, nz = rho - rc, z - zc
            size = np.maximum(np.hypot(nr, nz), 1e-12)
            nr, nz = nr / size, nz / size
        return radial * nr[:, None] + np.outer(nz, self.d)

    def point(self, rho, z, u):
        return self.a + z * self.d + rho * (math.cos(u) * self.e1 + math.sin(u) * self.e2)


class Sphere:
    def __init__(self, c, r):
        self.c, self.r, self.scale = np.asarray(c, float), float(r), float(r)

    kind = "sphere"

    def signed(self, p):
        return np.linalg.norm(np.atleast_2d(p) - self.c, axis=1) - self.r

    def normal(self, p):
        v = np.atleast_2d(p) - self.c
        return v / np.maximum(np.linalg.norm(v, axis=1), 1e-12)[:, None]


_mesh_tol = TOL


def _tol(model=None):
    return _mesh_tol


# ---------------------------------------------------------------- fitting

def _circle2d(x, y):
    if len(np.unique(np.round(np.c_[x, y], 5), axis=0)) < 3:
        return None
    sol = np.linalg.lstsq(np.c_[x, y, np.ones_like(x)], x * x + y * y, rcond=None)[0]
    cx, cy = sol[0] / 2, sol[1] / 2
    r2 = sol[2] + cx * cx + cy * cy
    return (cx, cy, math.sqrt(r2)) if r2 > 0 else None


def _fits(model, P):
    return model is not None and np.abs(model.signed(P)).max() <= _tol(model)


def _cylinder(d, P):
    d, e1, e2 = _frame(d)
    c = _circle2d(P @ e1, P @ e2)
    if c is None or c[2] > 1e4:
        return None
    return Revolved(c[0] * e1 + c[1] * e2, d, line=(c[2], 0))


def _cone(N, P):
    """Cone through facets whose normals N all make the same angle with the axis."""
    if len(np.unique(np.round(N, 4), axis=0)) < 3:
        return None
    Nc = N - N.mean(axis=0)
    d = np.linalg.eigh(Nc.T @ Nc)[1][:, 0]
    return _line_on_direction(d, P)


def _line_on_direction(d, P):
    d, e1, e2 = _frame(d)
    x, y, z = P @ e1, P @ e2, P @ d
    sol = np.linalg.lstsq(np.c_[2 * x, 2 * y, np.ones_like(x), z, z * z], x * x + y * y, rcond=None)[0]
    cx, cy = sol[0], sol[1]
    rho = np.hypot(x - cx, y - cy)
    if np.ptp(z) < 1e-6:
        return None
    k, c0 = np.polyfit(z, rho, 1)
    if abs(k) > math.tan(math.radians(80)):
        return None                     # nearly flat: leave it as a plane
    if _straight(k, z):
        k, c0 = 0, rho.mean()
    return Revolved(cx * e1 + cy * e2, d, line=(c0, k))


def _straight(k, z):
    """Is a line profile's slope k (over heights z) no slope at all? (A rounded-off mesh
    tilts a cylinder's fit by a hair; straightened, it still fits within tolerance.)"""
    return abs(k) < 1e-4 or abs(k) * np.ptp(z) / 2 <= _tol()


def _rings_refit(model, P):
    """Refit a cylinder/cone from its circular vertex rings: much more precise than normals."""
    rho, z, _ = model.local(P)
    levels = np.round(z / 0.01)
    centers, radii, heights, normals = [], [], [], []
    for lv in np.unique(levels):
        ring = P[levels == lv]
        if len(ring) < 3:
            continue
        mid = ring.mean(axis=0)
        n = np.linalg.svd(ring - mid)[2][2]
        _, e1, e2 = _frame(n)
        c = _circle2d((ring - mid) @ e1, (ring - mid) @ e2)
        if c is None:
            continue
        centers.append(mid + c[0] * e1 + c[1] * e2)
        radii.append(c[2])
        heights.append(lv)
        normals.append(n if n @ model.d > 0 else -n)
    if len(centers) < 2:
        return None
    centers = np.array(centers)
    lo, hi = np.argmin(heights), np.argmax(heights)
    d = centers[hi] - centers[lo]
    if np.linalg.norm(d) < 1e-6:
        return None
    d /= np.linalg.norm(d)
    if model.line[1] == 0:
        return _cylinder(d, P)
    return _line_on_direction(d, P)


def _sphere(P):
    if len(P) < 4 or np.linalg.svd(P - P.mean(axis=0), compute_uv=False)[2] < 1e-6:
        return None
    sol = np.linalg.lstsq(np.c_[2 * P, np.ones(len(P))], (P * P).sum(axis=1), rcond=None)[0]
    r2 = sol[3] + sol[:3] @ sol[:3]
    return Sphere(sol[:3], math.sqrt(r2)) if 0 < r2 < 1e8 else None


def _on_axis(axis, P):
    """Cone/cylinder (line profile) or torus/sphere (circle profile) around a known axis."""
    base = Revolved(axis[0], axis[1], line=(1.0, 0))
    rho, z, _ = base.local(P)
    if len(np.unique(np.round(np.c_[rho, z], 5), axis=0)) < 3:
        return None
    if np.ptp(z) > 1e-6:
        k, c0 = np.polyfit(z, rho, 1)
        if _straight(k, z):
            k, c0 = 0, rho.mean()
        line = Revolved(axis[0], axis[1], line=(c0, k))
        if abs(k) < math.tan(math.radians(80)) and _fits(line, P):
            return line
    c = _circle2d(rho, z)
    if c is None or c[2] > 1e4:
        return None
    rc = 0.0 if abs(c[0]) < TOL else c[0]
    return Revolved(axis[0], axis[1], circle=(rc, c[1], c[2]))


# ---------------------------------------------------------------- growing patches

def _normal_spread(model, pts):
    n = model.normal(pts)
    return math.acos(max(-1.0, min(1.0, float((n @ n.T).min()))))


class Region:
    """A set of smoothly connected facets, with fast "which facets lie on this surface" tests."""

    def __init__(self, mesh, facets):
        self.mesh = mesh
        self.facets = np.asarray(facets)
        pos = {f: i for i, f in enumerate(self.facets)}
        lists = [mesh.fverts[f] for f in self.facets]
        lens = np.array([len(l) for l in lists])
        self.vids, self.vinv = np.unique(np.concatenate(lists), return_inverse=True)
        self.starts = np.r_[0, np.cumsum(lens)[:-1]]
        self.fv = np.split(self.vinv, np.cumsum(lens)[:-1])
        self.P = mesh.pts[self.vids]
        self.n = mesh.fn[self.facets]
        self.c = mesh.fcent[self.facets]
        self.nbrs = [[pos[g] for g in mesh.nbrs[f] if g in pos] for f in self.facets]
        self.size = np.array([np.ptp(self.P[v], axis=0).max() for v in self.fv])
        # Needle-thin facets have unreliable plane normals even when their corners are exact.
        self.sliver = mesh.farea[self.facets] < 0.05 * self.size ** 2
        self.free = np.ones(len(self.facets), bool)
        self.explored = np.zeros(len(self.facets), bool)

    def points(self, idx):
        return self.P[np.unique(np.concatenate([self.fv[i] for i in idx]))]

    def _test(self, model, idx):
        """For facets idx: does each lie on the surface, which way does it face, surface normals."""
        lens = np.array([len(self.fv[i]) for i in idx])
        verts = np.concatenate([self.fv[i] for i in idx])
        worst = np.maximum.reduceat(np.abs(model.signed(self.P[verts])), np.r_[0, np.cumsum(lens)[:-1]])
        surface_n = model.normal(self.c[idx])
        facing = np.einsum("ij,ij->i", self.n[idx], surface_n)
        if _loose:
            # the surface normal's own spread over the facet is allowed on top
            spread = np.array([_normal_spread(model, self.P[self.fv[i]]) for i in idx])
            aligned = np.abs(facing) >= np.cos(np.minimum(np.radians(LOOSE_NORMAL_DEG) + spread,
                                                          np.radians(60)))
            gap = LOOSE_MAX_GAP
        else:
            aligned = (np.abs(facing) >= math.cos(math.radians(NORMAL_DEG))) | self.sliver[idx]
            gap = MAX_SAG
        ok = ((worst <= _tol(model)) & (np.abs(model.signed(self.c[idx])) <= gap)
              & aligned & self.free[idx])
        return ok, facing > 0, surface_n

    def grow(self, model, seeds):
        """Facets connected to the seeds that lie on the model surface, and whether it's convex.

        A facet only joins if it bends away from its neighbour the way the surface
        does. (A flat wall beside a rounded edge can have its corners on some huge
        circle, but it meets the rounded edge without the bend that circle predicts.)
        """
        seeds = list(seeds)
        ok, out, sn = self._test(model, seeds)
        if not ok.all() or (out != out[0]).any():
            return None, None
        convex = bool(out[0])
        slack = math.radians(BEND_SLACK_DEG)
        small = 0.1 * model.scale
        normal = dict(zip(seeds, sn))

        def bends_right(i, j, nj):
            if (self.size[i] < small and self.size[j] < small) or self.sliver[i] or self.sliver[j]:
                return True             # tiny or needle-thin facets: the corner check is what counts
            actual = math.acos(min(1.0, self.n[i] @ self.n[j]))
            predicted = math.acos(min(1.0, normal[i] @ nj))
            return abs(actual - predicted) <= slack

        for a in seeds[1:]:
            if not bends_right(seeds[0], a, normal[a]):
                return None, None
        # A facet that fails only the bend test is looked at again each time another
        # of its neighbours joins, since the verdict depends on which neighbour it's
        # compared with.
        tested, frontier = {}, seeds
        while frontier:
            cand = sorted({j for i in frontier for j in self.nbrs[i]
                           if j not in normal and tested.get(j, (True,))[0]})
            new = [j for j in cand if j not in tested]
            if new:
                tested.update(zip(new, zip(*self._test(model, new))))
            frontier = []
            for j in cand:
                good, o, nj = tested[j]
                if good and o == convex and any(i in normal and bends_right(i, j, nj)
                                                 for i in self.nbrs[j]):
                    normal[j] = nj
                    frontier.append(j)
        return np.array(sorted(normal)), convex


def _refit(model, P, axis_fixed):
    if isinstance(model, Sphere):
        return _sphere(P)
    if axis_fixed:
        return _on_axis((model.a, model.d), P)
    if model.line:
        return _rings_refit(model, P)
    return None


def _grow_refit(region, model, seeds, axis_fixed):
    idx, convex = region.grow(model, seeds)
    if idx is None:
        return None
    for _ in range(8):  # each refit uses more of the surface, so it can reach a little further
        better = _refit(model, region.points(idx), axis_fixed)
        if better is None or (not axis_fixed and better.kind != model.kind):
            break
        idx2, convex2 = region.grow(better, seeds)
        if idx2 is None or len(idx2) < len(idx):
            break
        grew = len(idx2) > len(idx)
        model, idx, convex = better, idx2, convex2
        if not grew:
            break
    return model, idx, convex


def _profile_param(model, P):
    """Position of points along the profile: height for a line, angle round the tube for a circle."""
    rho, z, _ = model.local(P)
    if model.line:
        return z, 1.0
    rc, zc, r = model.circle
    return np.arctan2(z - zc, rho - rc), r


def _trim(region, model, idx):
    """Drop facets that leaked past a patch's true end circles.

    Where two surfaces meet tangentially, the first row of tiny facets on the
    neighbour can sit within tolerance of this surface too. A real end circle
    carries many mesh corners; a leaked sliver adds only a few stray ones, so the
    ends are taken as the outermost well-populated levels.
    """
    tol = _tol(model)
    for _ in range(3):
        verts = [region.fv[i] for i in idx]
        t, scale = _profile_param(model, region.P)
        used = np.unique(np.concatenate(verts))
        levels, counts = np.unique(np.round(t[used] * scale / tol), return_counts=True)
        busy = levels[counts >= max(4, 0.2 * counts.max())] * tol / scale
        if len(busy) == 0:
            break
        lo, hi = busy.min() - 2 * tol / scale, busy.max() + 2 * tol / scale
        keep = np.array([((t[v] >= lo) & (t[v] <= hi)).all() for v in verts])
        if keep.all():
            break
        idx = idx[keep]
        if len(idx) < 3:
            break
    return idx


def _candidate(mesh, region, model, seeds, axis_fixed=False):
    """Grow a hypothesis into a patch; return (Feature, facet positions, area) or None.

    Facets of a patch that grew but was rejected go into region.explored, so the
    same surface isn't regrown (and rejected again) from each of its facets.
    """
    if model is None:
        return None
    grown = _grow_refit(region, model, seeds, axis_fixed)
    if grown is None:
        return None
    model, idx, convex = grown
    options = []
    if isinstance(model, Revolved) and not (model.circle and model.circle[0] == 0):
        cut = _trim(region, model, idx)
        if len(cut) >= 3 and len(cut) < len(idx):
            better = _refit(model, region.points(cut), axis_fixed)
            if (better is not None and (axis_fixed or better.kind == model.kind)
                    and region._test(better, cut)[0].all()):
                options.append((better, cut))
            else:
                options.append((model, cut))
    options.append((model, idx))
    result = None
    min_facets = 2 if _loose else 3
    for m, ids in options:
        if len(ids) < min_facets or len(np.unique(np.round(region.n[ids], 3), axis=0)) < min_facets:
            continue
        feature = _feature(mesh, m, region.facets[ids], convex)
        if feature is None:
            continue
        found = (feature, ids, float(mesh.farea[region.facets[ids]].sum()))
        if feature.kind != "trimmed":
            return found                # a clean outline: take it
        result = result or found
    if result is None and len(idx) >= 6:
        region.explored[idx] = True
    return result


def _best(candidates):
    """The candidate covering the most area; among those covering about as much, the one
    fitting its facets far more closely (two rings of points lie on a cone and on a
    sphere alike; only the true surface also runs close between them)."""
    found = [c for c in candidates if c]
    if not found:
        return None
    best = max(found, key=lambda c: c[2])
    close = [c for c in found if c[2] >= 0.99 * best[2]]
    tight = min(close, key=lambda c: c[0].worst)
    return tight if tight[0].worst < 0.5 * best[0].worst else best


def _seed_candidates(mesh, region, i, j):
    """Every surface the facet pair (i, j) could belong to, grown as far as it goes."""
    P = region.points([i, j])
    near = [i, j] + [k for k in region.nbrs[i] + region.nbrs[j] if region.free[k]]
    near_pts = region.points(near)
    found = []
    # The cylinder's axis: along the edge the two facets share (a strip cut straight
    # across), or square to both their normals (every facet of a cylinder faces straight
    # out from the axis, however the strip is cut into triangles).
    axes = [mesh.shared_direction(region.facets[i], region.facets[j])]
    cross = np.cross(region.n[i], region.n[j])
    if np.linalg.norm(cross) > math.sin(math.radians(1)):
        cross /= np.linalg.norm(cross)
        if axes[0] is None or abs(cross @ axes[0]) < math.cos(math.radians(0.5)):
            axes.append(cross)
    for d in axes:
        if d is None:
            continue
        cyl = _cylinder(d, P)
        if _fits(cyl, P):
            found.append(_candidate(mesh, region, cyl, [i, j]))
    for k in region.nbrs[j]:
        if k != i and region.free[k]:
            pts = region.points([i, j, k])
            cone = _cone(region.n[[i, j, k]], pts)
            if _fits(cone, pts):
                found.append(_candidate(mesh, region, cone, [i, j, k]))
            break
    sv = np.linalg.svd(region.n[near], compute_uv=False) if len(near) >= 3 else None
    if not any(found) and sv is not None and sv[2] <= 0.1 * sv[0]:
        # Irregular triangles (a round whose edges run at different heights) fool both
        # quick guesses: estimate the axis from the whole neighbourhood's normals and
        # let a least-squares fit through the corners settle it. (Only where the normals
        # all lie square to one direction, as on a cylinder; not on doubly curved areas.)
        axis = np.linalg.svd(region.n[near])[2][2]
        cyl = _cylinder(axis, near_pts)
        if cyl is not None:
            cyl = _tube_fit(near_pts, cyl)
            if _fits(cyl, near_pts) and cyl.line[0] < 1e3:
                found.append(_candidate(mesh, region, cyl, [i, j]))
    # The same axis may carry a torus or sphere that explains more than one ring of facets.
    for c in [c for c in found if c]:
        ring = _on_axis((c[0].model.a, c[0].model.d), near_pts)
        if _fits(ring, near_pts) and ring.circle:
            found.append(_candidate(mesh, region, ring, [i, j], axis_fixed=True))
    for pts in (near_pts, P):
        ball = _sphere(pts)
        if _fits(ball, pts):
            found.append(_candidate(mesh, region, ball, [i, j]))
            break
    return _best(found)


def _whole_region(mesh, region):
    """The one cylinder, cone or sphere a whole smooth region lies on, if it is all one
    surface (a countersink, a plain hole or boss between two creases). Fitted from every
    facet at once: a cone cut into long thin triangles fools any guess from a few facets
    (a thin strip of cone looks like a tilted cylinder)."""
    idx = np.arange(len(region.facets))
    P, N = region.points(idx), region.n
    if len(P) < MIN_CORNERS["cylinder"]:
        return None
    guesses = []
    sv, V = np.linalg.svd(N, full_matrices=False)[1:]
    if sv[2] <= 0.1 * sv[0]:
        guesses.append(_cylinder(V[2], P))      # every normal square to the axis
    guesses.append(_cone(N, P))                  # every normal at one angle to the axis
    found = []
    for model in guesses:
        if model is None:
            continue
        model = _tube_fit(P, model)
        if _fits(model, P):
            seed = [int(np.argmax(mesh.farea[region.facets]))]
            found.append(_candidate(mesh, region, model, seed))
        elif model is not None:
            found.append(_whole_despite_seams(mesh, region, model))
    ball = _sphere(P)
    if _fits(ball, P):
        found.append(_candidate(mesh, region, ball, [0]))
    return _best(found)


def _whole_despite_seams(mesh, region, model):
    """The whole region on this surface although a few corners are off it: meshes where
    one face's triangulation put extra points on the straight chords of another face's
    outline (a few microns inside the curve). Accepted only if nearly every corner lies
    on the surface and the stray ones are all on the region's outline, close by."""
    idx = np.arange(len(region.facets))
    P = region.P
    for _ in range(4):
        off = np.abs(model.signed(P))
        keep = off <= max(_tol(model), 3 * np.median(off))
        if keep.mean() < 0.5 or keep.sum() < MIN_CORNERS[model.kind]:
            return None
        model = _tube_fit(P[keep], model)
        if model is None:
            return None
    off = np.abs(model.signed(P))
    on = off <= _tol(model)
    if on.mean() < WHOLE_ON_SURFACE or off.max() > WHOLE_STRAY:
        return None
    T = mesh.tris[np.concatenate([mesh.ftris[f] for f in region.facets])]
    uniq, counts = np.unique(_edges(T), axis=0, return_counts=True)
    if not set(region.vids[~on].tolist()) <= set(uniq[counts == 1].ravel().tolist()):
        return None
    ok, out, surface_n = region._test(model, idx)
    facing = np.abs(np.einsum("ij,ij->i", region.n, surface_n))
    if (out != out[0]).any() or (facing < math.cos(math.radians(NORMAL_DEG))).any():
        return None
    feature = _feature(mesh, model, region.facets, bool(out[0]))
    if feature is None:
        return None
    return feature, idx, float(mesh.farea[region.facets].sum())


def find_features(pts, tris):
    """Return a Feature for every curved patch that can be rebuilt exactly."""
    return analyze(pts, tris)[1]


def analyze(pts, tris):
    """Return (mesh, features, tolerance): the facet structure and every rebuildable patch."""
    global _mesh_tol
    mesh = Mesh(pts, tris)
    _mesh_tol = TOL + FLOAT_TOL * float(np.abs(pts).max())
    regions = [Region(mesh, g) for g in mesh.smooth_regions()]
    features = []

    # Pass 0: smooth regions that are all one surface, fitted whole.
    for region in regions:
        if len(region.facets) >= 6:
            best = _whole_region(mesh, region)
            if best:
                feature, idx, _ = best
                region.free[idx] = False
                features.append(feature)

    # Screw threads, among the curved facets still unexplained: before any pass that
    # tries surfaces round known axes (a thread's crest, root and flanks lie close to a
    # cylinder or cone on its own axis, which is often a hole's too) or the loose pass
    # (which would take pieces of them for cylinders and cones)
    import threads
    free = [region.facets[region.free] for region in regions]
    screws = []
    for thread, fids in threads.find(mesh, np.concatenate(free) if free else np.zeros(0, int)):
        screws += threads.patches(mesh, thread, fids)
    claimed = [int(f) for x in screws for f in x.facets]
    for region in regions:
        region.free &= ~np.isin(region.facets, claimed)

    # Pass 0b: surfaces on the axes found so far (a lug round its screw hole), before
    # pairs of facets can propose surfaces of their own: on a coarse, rounded-off mesh a
    # strip of three or four facets fits a cylinder of almost any radius and tilt.
    features += _axis_pass(mesh, regions, _distinct_axes(features, mesh), AXIS_NOISE * mesh.noise)

    # Pass 1: from each pair of neighbouring facets, work out what surface they're on.
    for region in regions:
        for i in range(len(region.facets)):
            for j in region.nbrs[i][:3]:
                if not (region.free[i] and region.free[j]) or region.explored[i]:
                    continue
                best = _seed_candidates(mesh, region, i, j)
                if best:
                    feature, idx, _ = best
                    region.free[idx] = False
                    features.append(feature)
                    break

    # Pass 2: tori (and anything missed) around the axes of the patches found so far.
    # A curved rounded edge always shares its axis with a neighbouring hole, pin or
    # rounded corner, so each new patch can unlock its neighbours.
    done = []
    while True:
        # each round only tries the axes found since the last one
        axes = [x for x in _distinct_axes(features, mesh)
                if not any(abs(x[1] @ d) > math.cos(math.radians(0.5))
                           and np.linalg.norm(np.cross(x[0] - a, d)) < TOL * 4 for a, d, _ in done)]
        if not axes:
            break
        done += axes
        # (on a known axis only the radius is free, so a file that rounded its corners
        # can be allowed that much more: the CAD program's own export is often as loose)
        features += _axis_pass(mesh, regions, axes, AXIS_NOISE * mesh.noise)
    # Pass 3: fillets between flat faces that the passes above couldn't make out (cut
    # into few or irregular strips): a cylinder touching both faces, only its radius to find
    import fillets
    taken = np.zeros(len(mesh.farea), bool)
    for f in features:
        taken[f.facets] = True
    rolled = fillets.find(mesh, regions, taken)
    claimed = [int(f) for x in rolled for f in x.facets]
    for region in regions:
        region.free &= ~np.isin(region.facets, claimed)
    features += rolled
    # a thread usually ends in a countersink or chamfer round its own axis, which the
    # passes above couldn't try (the thread wasn't known yet). Held to the thread's own
    # tolerance: the CAD program drew them together, both a hundredth or so off true.
    axes = [(t.a, t.d, 2 * float(t.profile[:, 1].max()) + 3) for t in {id(p.model.thread): p.model.thread for p in screws}.values()]
    ends = _axis_pass(mesh, regions, axes, threads.ON_TOL + mesh.noise) if axes else []
    features += _loose_pass(mesh, regions, features)
    features = _band_tori(mesh, regions, features)
    # (the thread's pieces ahead of what was found round its ends: if the solid won't
    # check out with everything, the search for culprits keeps what comes first)
    return mesh, _merge_same_surface(mesh, features) + screws + ends, _mesh_tol


def _tube_fit(P, model, steps=30):
    """Least-squares torus, cylinder or cone through points P, starting from `model`.

    Gauss-Newton on the fewest parameters: the axis tilted by two small angles, the
    centre moved across it (and along it for a torus), and the radii (for a cone, the
    radius where the axis point is and how fast it grows).
    """
    torus = not model.line
    cone = not torus and model.line[1] != 0
    d0, e1, e2 = _frame(model.d)
    if torus:
        rc, zc, r = model.circle
        c0, x = model.a + zc * d0, np.array([0, 0, 0, 0, 0, rc, r], float)
    elif cone:
        c0, x = model.a, np.array([0, 0, 0, 0, model.line[0], model.line[1]], float)
    else:
        c0, x = model.a, np.array([0, 0, 0, 0, model.line[0]], float)

    def res(x):
        d = d0 + x[0] * e1 + x[1] * e2
        d /= np.linalg.norm(d)
        c = c0 + x[2] * e1 + x[3] * e2 + (x[4] * d0 if torus else 0)
        v = P - c
        z = v @ d
        rho = np.linalg.norm(v - np.outer(z, d), axis=1)
        if torus:
            return np.hypot(rho - x[5], z) - x[6]
        if cone:
            return (rho - x[4] - x[5] * z) / math.hypot(1, x[5])
        return rho - x[4]

    f = res(x)
    for _ in range(steps):
        J = np.empty((len(P), len(x)))
        for k in range(len(x)):
            h = 1e-7 * max(1.0, abs(x[k]))
            dx = x.copy()
            dx[k] += h
            J[:, k] = (res(dx) - f) / h
        step = np.linalg.lstsq(J, -f, rcond=None)[0]
        x2 = x + step
        f2 = res(x2)
        while f2 @ f2 > f @ f and np.abs(step).max() > 1e-14:
            step /= 2
            x2 = x + step
            f2 = res(x2)
        if f2 @ f2 > f @ f:
            break
        x, f = x2, f2
        if np.abs(step).max() < 1e-12:
            break
    if not np.all(np.isfinite(x)):
        return None
    d = d0 + x[0] * e1 + x[1] * e2
    d /= np.linalg.norm(d)
    c = c0 + x[2] * e1 + x[3] * e2
    if torus:
        c = c + x[4] * d0
        return Revolved(c, d, circle=(abs(x[5]), float(c @ d), abs(x[6])))
    if cone:
        # the model measures heights from the axis point nearest the origin
        return Revolved(c, d, line=(float(x[4] - x[5] * (c @ d)), float(x[5])))
    return Revolved(c, d, line=(abs(x[4]), 0))


def _band_tori(mesh, regions, features):
    """Replace chains of small pieces by the torus (or cylinder) they approximate.

    A rounded edge that follows a curve, with no neighbouring hole or pin to supply its
    axis, is easily cut into many small pieces of the fillet radius: short cylinder
    strips (the idea of the stlToSolid project) or, on coarse meshes, little "spheres"
    (a short stretch of a torus fits a sphere of its tube radius too). The pieces' centres
    run along the torus's centre circle, so consecutive pieces are fitted as one torus,
    and the run is regrown as a single patch.
    """
    global _loose

    def piece(f):
        """(radius, point on the tube's centre line) for a small piece, else None."""
        m = f.model
        if isinstance(m, Sphere):
            P = mesh.pts[np.unique(np.concatenate([mesh.fverts[x] for x in f.facets]))]
            return (m.r, m.c) if np.ptp(P, axis=0).max() <= 2.5 * m.r else None
        if isinstance(m, Revolved) and m.line and m.line[1] == 0 and f.span < TWO_PI:
            if f.hi - f.lo <= 4 * m.line[0]:
                return m.line[0], m.a + (f.lo + f.hi) / 2 * m.d
        return None

    small = {k: piece(f) for k, f in enumerate(features)}
    small = {k: v for k, v in small.items() if v is not None}
    if len(small) < 3:
        return features
    keys = list(small)
    verts = {k: set(np.concatenate([mesh.fverts[x] for x in features[k].facets]).tolist()) for k in keys}
    link = {k: [] for k in keys}
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            ra, rb = small[a][0], small[b][0]
            if abs(ra - rb) <= max(0.05 * ra, 2 * _tol()) and verts[a] & verts[b]:
                link[a].append(b)
                link[b].append(a)
    region_of = {}
    for region in regions:
        for pos, f in enumerate(region.facets):
            region_of[int(f)] = (region, pos)

    def ordered(chain):
        """The chain's pieces in order along it (a walk from one end)."""
        members = set(chain)
        ends = [k for k in chain if len([j for j in link[k] if j in members]) <= 1]
        k = ends[0] if ends else chain[0]
        out, seen = [k], {k}
        while True:
            nxt = [j for j in link[k] if j in members and j not in seen]
            if not nxt:
                break
            # the nearest neighbour keeps the walk on the chain where it branches
            k = min(nxt, key=lambda j: np.linalg.norm(small[j][1] - small[out[-1]][1]))
            out.append(k)
            seen.add(k)
        return out

    def points(ks):
        return mesh.pts[np.unique(np.concatenate([mesh.fverts[x] for k in ks for x in features[k].facets]))]

    def tube(ks):
        """The torus or cylinder through these pieces, if their mesh corners lie on it."""
        C = np.array([small[k][1] for k in ks])
        r = float(np.mean([small[k][0] for k in ks]))
        mid = C.mean(axis=0)
        vt = np.linalg.svd(C - mid)[2]
        guesses = [Revolved(mid, vt[0], line=(r, 0))]           # a straight fillet
        if len(ks) >= 3:
            _, e1, e2 = _frame(vt[2])
            c = _circle2d((C - mid) @ e1, (C - mid) @ e2)
            if c is not None and c[2] < 1e3 * r:
                centre = mid + c[0] * e1 + c[1] * e2
                guesses.append(Revolved(centre, vt[2], circle=(c[2], float(centre @ _frame(vt[2])[0]), r)))
        P = points(ks)
        fitted = [m for m in (_tube_fit(P, g) for g in guesses) if m is not None]
        fitted = [(np.abs(m.signed(P)).max(), k, m) for k, m in enumerate(fitted)]
        if not fitted:
            return None
        err, _, model = min(fitted)
        if err > 2 * _tol():
            return None
        if model.circle and model.circle[0] <= model.circle[2]:
            return None                                    # not a tube round a hole in the middle
        return model

    drop, added, seen = set(), [], set()
    for start in keys:
        if start in seen or len(link[start]) == 0:
            continue
        chain, stack = [], [start]
        seen.add(start)
        while stack:
            k = stack.pop()
            chain.append(k)
            for j in link[k]:
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
        if len(chain) < 3:
            continue
        run_order = ordered(chain)
        i = 0
        while i + 3 <= len(run_order):
            # the longest run from i that one tube explains
            j, model = i + 3, tube(run_order[i:i + 3])
            if model is None:
                i += 1
                continue
            while j < len(run_order):
                better = tube(run_order[i:j + 1])
                if better is None:
                    break
                model, j = better, j + 1
            run = run_order[i:j]
            i = j
            facets = np.concatenate([features[k].facets for k in run])
            if int(facets[0]) not in region_of:
                continue
            region, seed = region_of[int(facets[0])]
            positions = [region_of[int(f)][1] for f in facets
                         if int(f) in region_of and region_of[int(f)][0] is region]
            region.free[positions] = True
            _loose = True
            try:
                got = _candidate(mesh, region, model, [seed], axis_fixed=True)
            finally:
                _loose = False
            area = float(mesh.farea[facets].sum())
            if got and got[0].model.kind == model.kind and got[2] >= 0.9 * area:
                feature, idx, _ = got
                region.free[idx] = False
                drop.update(run)
                # pieces the new patch swallowed whole go too
                taken = set(region.facets[idx].tolist())
                drop.update(k for k in keys if set(features[k].facets.tolist()) <= taken)
                added.append(feature)
            else:
                region.free[positions] = False
    if not added:
        return features
    return [f for k, f in enumerate(features) if k not in drop] + added


def _loose_pass(mesh, regions, features):
    """Pass 3: whatever curved facets are left, accept any surface that explains them."""
    global _loose
    axes = _distinct_axes(features, mesh)
    found_all = []
    _loose = True
    try:
        for region in regions:
            region.explored[:] = False
            for i in range(len(region.facets)):
                for j in region.nbrs[i][:3]:
                    if not (region.free[i] and region.free[j]) or region.explored[i]:
                        continue
                    found = [_seed_candidates(mesh, region, i, j)]
                    near = [i, j] + [k for k in region.nbrs[i] + region.nbrs[j] if region.free[k]]
                    pts = region.points(near)
                    for axis in _near(axes, pts):
                        m = _on_axis(axis, pts)
                        if _fits(m, pts):
                            found.append(_candidate(mesh, region, m, [i, j], axis_fixed=True))
                    best = _best(found)
                    if best:
                        feature, idx, _ = best
                        region.free[idx] = False
                        found_all.append(feature)
                        break
    finally:
        _loose = False
    return found_all


def _simple_outline(mesh, fids):
    """True if no vertex appears more than once on the patch's outline loops."""
    T = np.concatenate([mesh.tris[mesh.ftris[f]] for f in fids])
    uniq, counts = np.unique(_edges(T), axis=0, return_counts=True)
    rim = uniq[counts == 1].ravel()
    return len(rim) == 0 or np.bincount(rim).max() <= 2


# ---------------------------------------------------------------- snapping

def _circle_fixed_r(x, y, r, cx, cy):
    """Centre of the radius-r circle best fitting points (x, y), starting from (cx, cy)."""
    for _ in range(30):
        dx, dy = x - cx, y - cy
        d = np.maximum(np.hypot(dx, dy), 1e-12)
        step = np.linalg.lstsq(np.c_[-dx / d, -dy / d], -(d - r), rcond=None)[0]
        cx, cy = cx + step[0], cy + step[1]
        if np.abs(step).max() < 1e-13:
            break
    return cx, cy


def _with_radius(model, P, r):
    """The same kind of surface with radius r, repositioned to fit points P best."""
    if isinstance(model, Sphere):
        c = model.c
        for _ in range(30):
            v = P - c
            d = np.maximum(np.linalg.norm(v, axis=1), 1e-12)
            step = np.linalg.lstsq(-v / d[:, None], -(d - r), rcond=None)[0]
            c = c + step
            if np.abs(step).max() < 1e-13:
                break
        return Sphere(c, r)
    if model.line and model.line[1] == 0:           # cylinder: keep the direction
        x, y = P @ model.e1, P @ model.e2
        cx, cy = _circle_fixed_r(x, y, r, model.a @ model.e1, model.a @ model.e2)
        return Revolved(cx * model.e1 + cy * model.e2, model.d, line=(r, 0))
    if model.circle:                                # torus or cap: keep the axis
        rho, z, _ = model.local(P)
        rc, zc, _ = model.circle
        if rc == 0:                                 # sphere cap: slide the centre along the axis
            for _ in range(30):
                d = np.maximum(np.hypot(rho, z - zc), 1e-12)
                J = (-(z - zc) / d)[:, None]
                step = float(np.linalg.lstsq(J, -(d - r), rcond=None)[0][0])
                zc += step
                if abs(step) < 1e-13:
                    break
            return Revolved(model.a, model.d, circle=(0.0, zc, r))
        rc, zc = _circle_fixed_r(rho, z, r, rc, zc)
        return Revolved(model.a, model.d, circle=(rc, zc, r))
    return None


def _radius(model):
    if isinstance(model, Sphere):
        return model.r
    if not isinstance(model, Revolved):
        return None
    if model.line:
        return model.line[0] if model.line[1] == 0 else None
    return model.circle[2]


def _rebuilt(mesh, feature, model):
    """The feature on a new surface, if that surface still explains its facets."""
    global _loose
    for loose in (False, True):
        _loose = loose
        try:
            f = _feature(mesh, model, feature.facets, feature.convex)
        finally:
            _loose = False
        if f is not None:
            return f
    return None


def half_rings(mesh, feature):
    """A patch that runs all the way round its axis but is cut to shape (its outline
    wraps the surface's seam) as two half rings on the same surface, or None."""
    model = feature.model
    if feature.kind != "trimmed" or not isinstance(model, Revolved) or feature.span < TWO_PI:
        return None
    centres = np.array([mesh.pts[np.unique(mesh.fverts[f])].mean(axis=0) for f in feature.facets])
    u = model.angle(model.local(centres)[2])
    global _loose, _anchored
    out = []
    for side in (np.cos(u) >= 0, np.cos(u) < 0):
        _loose, _anchored = True, True     # the whole ring has already passed the checks
        try:
            half = _feature(mesh, model, feature.facets[side], feature.convex) if side.any() else None
        finally:
            _loose, _anchored = False, False
        if half is None:
            return None
        out.append(half)
    return out


def snap(mesh, features, round_unit=None):
    """Design intent: equal radii made exactly equal, near-axis-aligned axes made exact,
    and (with round_unit "mm" or "inch") radii set to round values. Each change is kept
    only if the patch's mesh corners still lie on the surface; otherwise it is undone."""
    from sizing import roundness, INCH
    tol = 2 * _tol()
    points = [mesh.pts[np.unique(np.concatenate([mesh.fverts[f] for f in x.facets]))] for x in features]
    changed = 0

    def fits(model, P):
        return model is not None and np.abs(model.signed(P)).max() <= tol

    # axes within a hair of X, Y or Z become exactly X, Y or Z
    for k, f in enumerate(features):
        m = f.model
        if not isinstance(m, Revolved):
            continue
        axis = np.eye(3)[np.argmax(np.abs(m.d))] * np.sign(m.d[np.argmax(np.abs(m.d))])
        if abs(m.d @ axis) >= 1 - 1e-12 or abs(m.d @ axis) < math.cos(math.radians(0.1)):
            continue
        if m.line and m.line[1] == 0:
            new = _cylinder(axis, points[k])
            new = _with_radius(new, points[k], m.line[0]) if new else None
        else:
            new = Revolved(m.a, axis, line=m.line, circle=m.circle)
        if fits(new, points[k]):
            g = _rebuilt(mesh, f, new)
            if g:
                features[k], changed = g, changed + 1

    # radii: group equal ones, then optionally round the group's value
    radii = [(k, _radius(f.model)) for k, f in enumerate(features)]
    radii = sorted([(r, k) for k, r in radii if r], key=lambda x: x[0])
    groups, current = [], []
    for r, k in radii:
        if current and r - current[-1][0] > 2e-3 * r:
            groups.append(current)
            current = []
        current.append((r, k))
    if current:
        groups.append(current)
    for group in groups:
        weights = np.array([mesh.farea[features[k].facets].sum() for _, k in group])
        target = float(np.average([r for r, _ in group], weights=weights))
        if round_unit:
            scale = INCH if round_unit == "inch" else 1.0
            for step in ((1 / 16, 0.05, 1 / 64, 0.01) if round_unit == "inch" else (0.5, 0.1)):
                nice = round(target / scale / step) * step * scale
                if nice > 0 and roundness(nice, round_unit) and abs(nice - target) <= 2e-3 * target:
                    target = nice
                    break
        for r, k in group:
            if r == target:
                continue
            new = _with_radius(features[k].model, points[k], target)
            if fits(new, points[k]):
                g = _rebuilt(mesh, features[k], new)
                if g:
                    features[k], changed = g, changed + 1
    return changed


def _merge_same_surface(mesh, features):
    """Join neighbouring patches that lie on the same surface (a hole split in two by a
    seam in the mesh, say), so each becomes one face."""
    verts = [set(np.concatenate([mesh.fverts[f] for f in x.facets]).tolist()) for x in features]
    merged = True
    while merged:
        merged = False
        for i in range(len(features)):
            for j in range(i + 1, len(features)):
                a, b = features[i], features[j]
                if a.convex != b.convex or not (verts[i] & verts[j]):
                    continue
                if (np.abs(a.model.signed(mesh.pts[list(verts[j])])).max() > _tol()
                        or np.abs(b.model.signed(mesh.pts[list(verts[i])])).max() > _tol()):
                    continue
                fids = np.concatenate([a.facets, b.facets])
                if not _simple_outline(mesh, fids):
                    continue                # outline touches itself (e.g. equal holes crossing)
                both = _feature(mesh, a.model, fids, a.convex)
                if both is None:
                    continue
                features[i] = both
                verts[i] |= verts.pop(j)
                features.pop(j)
                merged = True
                break
            if merged:
                break
    return features


def _distinct_axes(features, mesh=None):
    """[(point, direction, reach)]: each distinct axis, and how far from it a neighbouring
    rounded edge could plausibly lie (a few times the largest radius on that axis).

    A dome or ball contributes an axis through its centre square to each flat face it
    borders: a rounded edge between the two is a torus on that axis."""
    found = []
    for f in features:
        m = f.model
        if isinstance(m, Revolved):
            size = m.line[0] + abs(m.line[1]) * 10 if m.line else m.circle[0] + m.circle[2]
            found.append((m.a, m.d, size))
        elif isinstance(m, Sphere) and mesh is not None:
            inside = set(f.facets.tolist())
            P = mesh.pts[np.unique(np.concatenate([mesh.fverts[x] for x in f.facets]))]
            if np.ptp(P, axis=0).max() < m.r:
                continue                      # a small piece of something else
            border = {g for x in f.facets for g in mesh.nbrs[x] if g not in inside}
            big = [g for g in border if mesh.farea[g] > 0.05 * m.r * m.r]
            for g in big:
                found.append((m.c, mesh.fn[g], m.r))
    axes = []
    for a, d, size in found:
        a = a - (a @ d) * d
        for k, (b, e, reach) in enumerate(axes):
            if abs(d @ e) > math.cos(math.radians(0.5)) and np.linalg.norm(np.cross(a - b, e)) < TOL * 4:
                axes[k] = (b, e, max(reach, 2 * size + 3))
                break
        else:
            axes.append((a, d, 2 * size + 3))
    return axes


def _axis_pass(mesh, regions, axes, tol=None):
    """Patches on the given axes (cylinders, cones, tori) among the free facets, held to
    `tol` if given (else the mesh's own tolerance)."""
    global _mesh_tol, _axis_vouched
    strict = _mesh_tol
    if not tol or tol <= strict:
        return _axis_patches(mesh, regions, axes)
    # the wider allowance only for straight profiles (cylinders, cones: only the
    # radius and slope free); a torus's three free numbers find false fits in it
    _mesh_tol, _axis_vouched = tol, True
    try:
        out = _axis_patches(mesh, regions, axes, lines_only=True)
    finally:
        _mesh_tol, _axis_vouched = strict, False
    return out + _axis_patches(mesh, regions, axes)


def _axis_patches(mesh, regions, axes, lines_only=False):
    out = []
    for region in regions:
        region.explored[:] = False
        for i in range(len(region.facets)):
            near = [i] + [k for k in region.nbrs[i] if region.free[k]]
            if not region.free[i] or region.explored[i] or len(near) < 2:
                continue
            pts = region.points(near)
            found = []
            for axis in _near(axes, pts):
                model = _on_axis(axis, pts)
                if lines_only and model is not None and not model.line:
                    continue
                if _fits(model, pts):
                    found.append(_candidate(mesh, region, model, [i], axis_fixed=True))
            if _axis_vouched:
                # (the wider allowance lets a cone a degree off a cylinder fit a strip of
                # it; a real cone, a chamfer or countersink, tapers far more)
                found = [c for c in found if not c or not c[0].model.line
                         or c[0].model.line[1] == 0 or abs(c[0].model.line[1]) >= math.tan(math.radians(MIN_VOUCHED_CONE_DEG))]
            best = _best(found)
            if best:
                feature, idx, _ = best
                region.free[idx] = False
                out.append(feature)
    return out


def _near(axes, pts, most=6):
    """Only the axes that pass close enough to these points to matter (the nearest few)."""
    if not axes:
        return []
    key = id(axes)
    if _near.cache.get("key") != key:
        _near.cache = {"key": key, "A": np.array([a for a, _, _ in axes]),
                       "D": np.array([d for _, d, _ in axes]), "R": np.array([r for _, _, r in axes]),
                       "axes": axes}
    c = _near.cache
    dist = np.linalg.norm(np.cross(pts.mean(axis=0) - c["A"], c["D"]), axis=1)
    ok = np.nonzero(dist <= c["R"])[0]
    ok = ok[np.argsort(dist[ok])][:most]
    return [(c["A"][k], c["D"][k]) for k in ok]


_near.cache = {}


# ---------------------------------------------------------------- features

def _arc(angles, min_points=6):
    """(start, span) of the arc covered by a set of angles; span 2*pi for a full ring."""
    a = np.sort(np.mod(angles, TWO_PI))
    if len(a) < 2:
        return None
    a = a[np.r_[True, np.diff(a) > 1e-6]]   # merge repeats (the same angle on several rings)
    if len(a) < 2:
        return None
    gaps = np.diff(np.r_[a, a[0] + TWO_PI])
    g = gaps.argmax()
    if len(a) >= min_points and gaps[g] <= math.radians(SMOOTH_DEG + 5):
        return 0.0, TWO_PI
    return a[(g + 1) % len(a)], TWO_PI - gaps[g]


def _angle_gap(a, b):
    return np.abs(np.mod(a - b + math.pi, TWO_PI) - math.pi)


@dataclass
class Feature:
    model: object
    label: str          # what it is, in plain words
    detail: str
    convex: bool
    kind: str           # "revolve" (cylinder/cone/torus/sphere cap), "wedge" (sphere corner),
                        # "ball", "trimmed" (any of these surfaces, cut to an arbitrary outline)
                        # or "blend" (a smooth freeform face where no simple surface fits)
    facets: np.ndarray  # mesh facets this patch replaces
    change: float       # volume the exact surface adds compared with the facets (mm^3)
    tolerance: float    # how far the actual change may differ from that
    worst: float        # largest gap between a facet and the true surface (mm)
    u0: float = 0.0     # "revolve": angular range around the axis
    span: float = TWO_PI
    lo: float = 0.0     # "revolve": profile range (height for a line, tube angle for a circle)
    hi: float = 0.0
    planes: tuple = ()  # "wedge": planes through the centre bounding the corner
    parts: tuple = ()   # "blend": the pieces it replaced (put back if it can't be built)
    depth: int = 0      # "blend": how many times it has been cut in two to make it fit

    def describe(self):
        return f"{self.label:<22} {self.detail}"


def _feature(mesh, model, fids, convex):
    tids = np.concatenate([mesh.ftris[f] for f in fids])
    T = mesh.tris[tids]
    P = mesh.pts[np.unique(T)]
    uniq, counts = np.unique(_edges(T), axis=0, return_counts=True)
    rim_edges = uniq[counts == 1]
    rim = mesh.pts[np.unique(rim_edges)]

    # How far each facet sits from the true surface, sampled at its edge midpoints
    # (exact for the near-parabolic gap between a flat facet and a curved surface).
    corners = mesh.pts[T]
    s = model.signed(((corners + corners[:, [1, 2, 0]]) / 2).reshape(-1, 3)).reshape(-1, 3)
    area = mesh.tarea[tids]
    change = (-1 if convex else 1) * float((area * s.mean(axis=1)).sum())
    gain = float((area * np.abs(s).mean(axis=1)).sum())
    worst = float(np.abs(s).max())
    if worst < 1e-5:
        return None
    if _loose and not _anchored and (worst > LOOSE_MAX_GAP or not _loose_support(mesh, fids, model.kind)):
        return None
    base = dict(convex=convex, change=change, tolerance=0.5 * gain + 1e-3, facets=np.asarray(fids),
                worst=worst)

    if isinstance(model, Sphere):
        if len(rim) == 0:
            label = "ball" if convex else "spherical cavity"
            return Feature(model, label, f"dia {2 * model.r:.3f}", kind="ball", **base)
        mid = rim.mean(axis=0)
        n = np.linalg.svd(rim - mid)[2][2]
        if np.abs((rim - mid) @ n).max() <= _tol(model):
            # Flat circular edge: a dome, dimple or hemispherical tip, handled as a revolved cap.
            if ((P - mid) @ n).mean() < 0:
                n = -n
            cap = Revolved(model.c, n, circle=(0.0, float(model.c @ n), model.r))
            return _revolved_feature(mesh, cap, P, rim, rim_edges, base)
        directions = (P - model.c) / model.r
        spread = math.acos(max(-1.0, min(1.0, (directions @ directions.mean(axis=0)).min()
                                         / np.linalg.norm(directions.mean(axis=0)))))
        if not _loose and (len(P) < MIN_CORNERS["sphere"]
                           or 2 * spread < math.radians(MIN_SPAN_DEG)):
            return None
        planes = _great_circles(model, P, mesh.pts, rim_edges)
        if planes is None or len(planes) > 4:
            if not _loose and not _enough_for_trimmed(mesh, fids, P, rim_edges):
                return None
            label = "rounded area (sphere)" if convex else "spherical hollow"
            return Feature(model, label, f"sphere r {model.r:.3f}, cut to shape", kind="trimmed", **base)
        return Feature(model, "rounded corner", f"sphere r {model.r:.3f}", kind="wedge",
                       planes=tuple(planes), **base)
    return _revolved_feature(mesh, model, P, rim, rim_edges, base)


def _loose_support(mesh, fids, kind):
    """stlToSolid-style minimum evidence: enough facets, turning through enough angle."""
    n = mesh.fn[fids]
    if len(np.unique(np.round(n, 3), axis=0)) < LOOSE_MIN_FACETS[kind]:
        return False
    return (n @ n.T).min() <= math.cos(math.radians(LOOSE_MIN_TURN_DEG[kind]))


def _enough_for_trimmed(mesh, fids, P, free_edges):
    """Is a patch with an outline off its natural edges believable?

    Where it leaves those edges it must meet something else at a real crease (a hole
    breaking out through a sloped face, a knob cut by another feature). A smooth
    hand-over there means it's more likely a slice of some larger curved surface
    (one ring of a torus, say, which lies exactly on a sphere).
    """
    if len(P) < TRIMMED_MIN_CORNERS or len(fids) < TRIMMED_MIN_FACETS:
        return False
    n = mesh.fn[fids]
    if (n @ n.T).min() > math.cos(math.radians(TRIMMED_MIN_TURN_DEG)):
        return False
    # (Mostly, anyway: two equal holes crossing are tangent at a couple of points.)
    inside = set(np.asarray(fids).tolist())
    crease = math.cos(math.radians(CREASE_DEG))
    smooth = total = 0.0
    for a, b in free_edges:
        sides = mesh.edge_facets.get((min(a, b), max(a, b)), [])
        mine = [f for f in sides if f in inside]
        other = [f for f in sides if f not in inside]
        length = float(np.linalg.norm(mesh.pts[a] - mesh.pts[b]))
        total += length
        if mine and other and mesh.fn[mine[0]] @ mesh.fn[other[0]] > crease:
            smooth += length
    return smooth <= 0.25 * total


def _great_circles(model, P, pts, rim_edges):
    """Planes through the sphere centre that together make up a patch's outline."""
    tol = _tol(model)
    normals = []
    for i, j in rim_edges:
        m = np.cross(pts[i] - model.c, pts[j] - model.c)
        if np.linalg.norm(m) < 1e-9:
            continue
        m /= np.linalg.norm(m)
        if not any(abs(m @ q) > math.cos(math.radians(1)) for q in normals):
            normals.append(m)
    if not 1 <= len(normals) <= 8:
        return None
    # both ends of every outline edge must lie on the same plane
    on = np.abs((pts[rim_edges] - model.c) @ np.array(normals).T) <= tol
    if not (on[:, 0] & on[:, 1]).any(axis=1).all():
        return None
    planes = []
    for m in normals:
        side = (P - model.c) @ m
        if side.min() >= -tol:
            planes.append(m)
        elif side.max() <= tol:
            planes.append(-m)
        else:
            return None                 # patch is on both sides: not a simple corner
    return planes


def _revolved_feature(mesh, model, P, rim, rim_edges, base):
    tol = _tol(model)
    rho, z, w = model.local(P)
    off_axis = rho > tol
    around = _arc(model.angle(w[off_axis]))
    if around is None:
        return None
    u0, span = around
    r_rho, r_z, r_w = model.local(rim) if len(rim) else (np.zeros(0),) * 3
    r_u = model.angle(r_w) if len(rim) else np.zeros(0)

    if model.line:
        lo, hi = z.min(), z.max()
        at = lambda end: np.abs(r_z - end) <= tol
        c0, k = model.line
        if k != 0 and span >= TWO_PI:
            # A pointed cone tip can be closed off by one flat facet with no corner at
            # the tip itself; if the outline doesn't reach that end, run to the tip.
            apex = -c0 / k
            if abs(apex - lo) < abs(apex - hi) and not at(lo).any():
                lo = apex
            elif abs(apex - hi) < abs(apex - lo) and not at(hi).any():
                hi = apex
    else:
        rc, zc, r = model.circle
        ring = _arc(np.arctan2(z - zc, rho - rc))
        if ring is None or ring[1] >= TWO_PI - 1e-9:
            return None                 # a complete tube (O-ring): not handled
        lo, hi = ring[0], ring[0] + ring[1]
        v = np.arctan2(r_z - zc, r_rho - rc)
        at = lambda end: _angle_gap(v, end) * r <= tol
        if rc == 0 and span >= TWO_PI:
            # Same for the pole of a dome or dimple closed off by one flat facet.
            lo, hi = math.remainder(lo, TWO_PI), math.remainder(lo, TWO_PI) + (hi - lo)
            if not at(hi).any() and hi > 0:
                hi = math.pi / 2
            if not at(lo).any() and lo < 0:
                lo = -math.pi / 2

    # Every outline edge must run along one of the patch's natural boundary lines
    # (both of its ends on the same line), so the outline traces the patch exactly.
    lines = [at(lo), at(hi)]
    if span < TWO_PI:
        on_axis = r_rho <= tol
        lines += [(_angle_gap(r_u, u0) * r_rho <= tol) | on_axis,
                  (_angle_gap(r_u, u0 + span) * r_rho <= tol) | on_axis]
    trimmed = False
    if len(rim_edges):
        index = {v: i for i, v in enumerate(np.unique(rim_edges))}
        a = np.array([index[i] for i in rim_edges[:, 0]])
        b = np.array([index[j] for j in rim_edges[:, 1]])
        along = np.any([line[a] & line[b] for line in lines], axis=0)
        if not along.all():
            # the outline leaves the natural edges: only OK as a well-supported trimmed patch
            vouched = _axis_vouched and len(base["facets"]) >= TRIMMED_MIN_FACETS
            if not _loose and not vouched and not _enough_for_trimmed(mesh, base["facets"], P, rim_edges[~along]):
                return None
            trimmed = True
    if not _loose:
        if len(P) < MIN_CORNERS[model.kind]:
            return None                 # too few points to be sure of the shape
        min_span = math.radians(MIN_SPAN_DEG)
        if span < min_span or (not model.line and hi - lo < min_span):
            return None
        if model.kind == "cone" and not _creased(mesh, base["facets"], z, lo, hi, model, tol):
            return None

    if model.circle and any(model.circle[0] + model.circle[2] * math.cos(v) < -tol
                            for v in np.linspace(lo, hi, 9)):
        return None                     # profile would cross the axis
    label, detail = _describe(model, base["convex"], span, lo, hi)
    if trimmed:
        detail += ", cut to shape"
    return Feature(model, label, detail, kind="trimmed" if trimmed else "revolve",
                   u0=u0, span=span, lo=lo, hi=hi, **base)


def _creased(mesh, fids, z, lo, hi, model, tol):
    """Does a one-ring cone meet a neighbour at a real crease?

    A single ring of facets on any rounded surface is exactly a slice of a cone, so
    a one-ring cone is only trusted if at least one of its end circles is a crease
    (or an open edge) rather than a smooth continuation.
    """
    if ((np.abs(z - lo) > tol) & (np.abs(z - hi) > tol)).any():
        return True                     # several rings of facets: a genuine cone
    inside = set(fids.tolist())
    soft = set()
    for f in fids:
        for g in mesh.nbrs[f]:
            if g in inside:
                continue
            if math.degrees(math.acos(min(1.0, mesh.fn[f] @ mesh.fn[g]))) >= CREASE_DEG:
                continue
            for a, b in mesh.shared[(min(f, g), max(f, g))]:
                ez = model.local(mesh.pts[[a, b]])[1]
                for end in (lo, hi):
                    if (np.abs(ez - end) <= tol).all():
                        soft.add(end)
    c0, k = model.line
    ends = {end for end in (lo, hi) if c0 + k * end > tol}   # a pointed tip has no edge
    return not ends <= soft


def _describe(model, convex, span, lo, hi):
    full = span >= TWO_PI
    kind = model.kind
    if kind == "cylinder":
        r = model.line[0]
        if full:
            return ("pin" if convex else "hole"), f"dia {2 * r:.3f}, length {hi - lo:.2f}"
        return "rounded edge", f"r {r:.3f}, {math.degrees(span):.0f} deg, length {hi - lo:.2f}"
    if kind == "cone":
        c0, k = model.line
        angle = math.degrees(math.atan(abs(k)))
        dias = f"dia {2 * (c0 + k * lo):.3f} -> {2 * (c0 + k * hi):.3f}"
        if full:
            return ("cone / chamfered pin" if convex else "countersink / chamfer"), f"{dias}, {angle:.1f} deg"
        return "chamfer (cone)", f"{dias}, {angle:.1f} deg, {math.degrees(span):.0f} deg around"
    if kind == "torus":
        rc, zc, r = model.circle
        return "rounded edge (curved)", f"r {r:.3f} around a {2 * rc:.3f} circle, {math.degrees(span):.0f} deg"
    r = model.circle[2]
    height = r * (math.sin(hi) - math.sin(lo))
    return ("dome" if convex else "dimple"), f"sphere dia {2 * r:.3f}, height {height:.2f}"


def summarize(features):
    counts = Counter(f.label for f in features)
    return ", ".join(f"{label} x{n}" for label, n in counts.most_common())
