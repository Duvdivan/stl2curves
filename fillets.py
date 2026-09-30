"""
Fillets between flat faces, found the way a designer sees them.

Two flat faces that would meet at a sharp edge, rounded off with a fillet: the fillet
is a cylinder touching both faces (a ball of the fillet radius rolled along the
edge). Its axis runs parallel to the edge where the faces would meet, so given the two
faces the only unknown is the radius, and every mesh corner on the rounded strip pins
it down. This finds fillets however coarsely or irregularly they were cut into
triangles, even when neighbouring strips meet at a sharp-looking angle, which is where
fitting surfaces to the strips alone struggles. (Constant-radius "rolling ball" blend
recovery between primary surfaces, as in reverse-engineering work by Várady, Benkő,
Kós and Martin.)
"""
import math

import numpy as np

import features as F
from features import Revolved

MIN_FACE_AREA = 0.3     # mm^2: flat facets at least this big can be the faces a fillet joins
MIN_ANGLE_DEG = 30      # the faces must meet at least this far from flat...
MAX_ANGLE_DEG = 170     # ...and not be (nearly) parallel
MAX_GAP = 8             # mm: faces further apart than this aren't joined by one fillet
FACE_TO_STRIP = 3       # each face a fillet joins is at least this many times its biggest strip


def _adjacency(mesh):
    """Every facet's neighbours across any shared edge, however sharp the bend."""
    nbrs = [set() for _ in range(len(mesh.farea))]
    for fs in mesh.edge_facets.values():
        for a in fs:
            for b in fs:
                if a != b:
                    nbrs[a].add(b)
    return nbrs


def find(mesh, regions, taken=None):
    """Fillet features between pairs of flat faces. Marks their facets as used in the
    regions so later passes leave them alone."""
    nbrs = _adjacency(mesh)
    nf = len(mesh.farea)
    used = np.zeros(nf, bool) if taken is None else taken.copy()
    tol = 2 * F._tol()
    faces = [f for f in np.argsort(-mesh.farea) if mesh.farea[f] >= MIN_FACE_AREA and not used[f]]
    plane_pt = {f: mesh.pts[mesh.fverts[f][0]] for f in faces}
    verts = mesh.fverts

    def on_plane(f, P):
        return np.abs((P - plane_pt[f]) @ mesh.fn[f]) <= tol

    # facet bounding boxes, to find faces near one another
    lo = np.array([mesh.pts[v].min(axis=0) for v in verts])
    hi = np.array([mesh.pts[v].max(axis=0) for v in verts])

    def gap(x, y):
        return float(np.linalg.norm(np.maximum(0, np.maximum(lo[x] - hi[y], lo[y] - hi[x]))))

    # Every pair of flat faces close enough to be joined by a fillet and meeting at a
    # real angle; the biggest pairs first (real faces before the strips of a fillet).
    fa = np.array(faces)
    pairs = []
    for i, A in enumerate(faces):
        near = fa[i + 1:][np.all((lo[fa[i + 1:]] <= hi[A] + MAX_GAP) & (hi[fa[i + 1:]] >= lo[A] - MAX_GAP), axis=1)]
        for B in near:
            c = float(mesh.fn[A] @ mesh.fn[B])
            if math.cos(math.radians(MAX_ANGLE_DEG)) <= c <= math.cos(math.radians(MIN_ANGLE_DEG)):
                pairs.append((min(mesh.farea[A], mesh.farea[B]), A, int(B)))
    pairs.sort(key=lambda x: -x[0])

    found = []
    for _, A, B in pairs:
        if used[A] or used[B]:
            continue
        for X, Y in ((A, B), (B, A)):
            # the strip starts at facets touching one face and lying near the other
            strip = [g for g in nbrs[X] if g != Y and not used[g]
                     and gap(g, Y) <= MAX_GAP]
            # (a fillet's strips are clearly smaller than both faces, as _grow demands in
            # the end; the long facets of a finely cut tube are big, but so are all their
            # neighbours: skip those pairs without the fitting)
            small = min(mesh.farea[A], mesh.farea[B]) / FACE_TO_STRIP
            if not any(mesh.farea[g] <= small for g in strip):
                continue
            got = _fillet(mesh, nbrs, used, X, Y, strip, plane_pt, on_plane, tol) if strip else None
            if got is not None:
                feature, facets = got
                used[facets] = True
                found.append(feature)
                break
    return found


def _fillet(mesh, nbrs, used, A, B, strip, plane_pt, on_plane, tol):
    """The fillet joining faces A and B, grown from the strip facets next to A, if any:
    (feature, facets) or None."""
    n1, n2 = mesh.fn[A], mesh.fn[B]
    d = np.cross(n1, n2)
    d /= np.linalg.norm(d)
    # a point on the line where the two faces would meet
    p = np.linalg.solve(np.array([n1, n2, d]), [n1 @ plane_pt[A], n2 @ plane_pt[B], 0.0])
    counts = [len(mesh.fverts[g]) for g in strip]
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    V = mesh.pts[np.concatenate([mesh.fverts[g] for g in strip])]
    hA, hB = (V - plane_pt[A]) @ n1, (V - plane_pt[B]) @ n2
    Q = V - p
    Q -= np.outer(Q @ d, d)          # position across the edge line
    for sign in (1.0, -1.0):         # convex edge (ball inside the part) or concave (outside)
        inside = (sign * hA <= tol) & (sign * hB <= tol)
        inside_facet = np.minimum.reduceat(inside.astype(int), starts).astype(bool)
        pick = np.repeat(inside_facet, counts) & (np.abs(hA) > tol)
        if not pick.any():
            continue
        # The ball's centre moves along w as the radius grows (1 mm off each face per mm);
        # a corner v lies on the ball of radius r when |q - r w| = r: a quadratic in r.
        w = -sign * (n1 + n2) / (1 + n1 @ n2)
        q = Q[pick]
        qa, qb, qc = w @ w - 1.0, -2.0 * (q @ w), np.einsum("ij,ij->i", q, q)
        disc = qb * qb - 4 * qa * qc
        ok = disc >= 0
        if abs(qa) < 1e-12:
            radii = -qc[qb != 0] / qb[qb != 0]
        else:
            root = np.sqrt(np.maximum(disc[ok], 0))
            radii = np.r_[(-qb[ok] - root) / (2 * qa), (-qb[ok] + root) / (2 * qa)]
        radii = np.sort(radii[radii > 10 * tol])
        if not len(radii):
            continue
        # the most popular radius (the corners of one fillet all agree)
        groups = np.split(radii, np.nonzero(np.diff(radii) > max(tol, 1e-4 * radii[-1]))[0] + 1)
        groups.sort(key=len, reverse=True)
        for g_r in groups[:3]:
            r = float(np.median(g_r))
            err = np.abs(np.linalg.norm(Q - r * w, axis=1) - r)
            fit = np.maximum.reduceat(err, starts) <= tol
            seeds = [g for g, f in zip(strip, fit & inside_facet) if f]
            got = _grow(mesh, nbrs, used, A, B, p + r * w, d, r, sign, on_plane, tol, seeds) if seeds else None
            if got is not None:
                return got
    return None


def _grow(mesh, nbrs, used, A, B, a, d, r, sign, on_plane, tol, seeds):
    def fits(g):
        if used[g] or g == A or g == B:
            return False
        q = mesh.pts[mesh.fverts[g]] - a
        q -= np.outer(q @ d, d)
        if np.abs(np.linalg.norm(q, axis=1) - r).max() > tol:
            return False
        # it must face out of the ball for a convex fillet, into it for a concave one
        c = mesh.fcent[g] - a
        c -= (c @ d) * d
        return (mesh.fn[g] @ c) / max(np.linalg.norm(c), 1e-12) * sign > 0.5

    seeds = [g for g in seeds if fits(g)]
    if not seeds:
        return None
    group, stack = set(seeds), list(seeds)
    while stack:
        g = stack.pop()
        for h in nbrs[g]:
            if h not in group and fits(h):
                group.add(h)
                stack.append(h)
    facets = np.array(sorted(group))
    # The faces it joins must be real faces, clearly bigger than its strips (two facets of
    # a finely cut pin, with the facets between them, can fit some other cylinder too)
    if min(mesh.farea[A], mesh.farea[B]) < FACE_TO_STRIP * mesh.farea[facets].max():
        return None
    # it must actually join both faces, and bend in between (a single flat strip with its
    # edges on both faces is as likely a chamfer)
    if not any(B in nbrs[g] for g in facets) or not any(A in nbrs[g] for g in facets):
        return None
    P = mesh.pts[np.unique(np.concatenate([mesh.fverts[g] for g in facets]))]
    if (~(on_plane(A, P) | on_plane(B, P))).sum() < 2:
        return None
    model = Revolved(a, d, line=(r, 0.0))
    F._loose, F._anchored = True, True
    try:
        feature = F._feature(mesh, model, facets, sign > 0)
    finally:
        F._loose, F._anchored = False, False
    if feature is None:
        return None
    return feature, facets
