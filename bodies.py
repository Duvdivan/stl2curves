"""Separate bodies that touch along an edge.

Some STLs are several closed bodies that meet exactly along edges (a model exported as
separate pieces). Where they meet, four triangles share one edge. Sewing welds such
bodies into one tangled shell, so they are split apart here and converted one group at
a time, then fused.
"""

import numpy as np


def _union_find(n):
    parent = np.arange(n)

    def find(a):
        root = a
        while parent[root] != root:
            root = parent[root]
        while parent[a] != root:
            parent[a], a = root, parent[a]
        return root

    return parent, find


def _partners(pts, tris, t_list, u, v):
    """Pair up the triangles around edge (u, v): each with the first triangle reached by
    turning about the edge into its own material (against its normal)."""
    d = pts[v] - pts[u]
    d /= np.linalg.norm(d)
    ang, turn, ok = [], [], []
    for t in t_list:
        a, b, c = pts[tris[t]]
        n = np.cross(b - a, c - a)
        third = pts[[w for w in tris[t] if w != u and w != v][0]] - pts[u]
        w = third - (third @ d) * d                 # direction from the edge into the triangle
        if np.linalg.norm(n) < 1e-12 or np.linalg.norm(w) < 1e-9:
            n, w = np.array([1.0, 0, 0]), np.cross(d, [1.0, 0, 0]) + 1e-3   # flat sliver: any side
            ok.append(False)
        else:
            ok.append(True)
        n /= np.linalg.norm(n)
        w /= np.linalg.norm(w)
        e2 = np.cross(d, w)
        ang.append(w)
        turn.append(1.0 if (-n) @ e2 > 0 else -1.0)   # which way is "into the material"
    ref = ang[0]
    ref2 = np.cross(d, ref)
    theta = [np.arctan2(w @ ref2, w @ ref) for w in ang]
    # which way each triangle runs along the edge: a body's two faces run opposite ways
    forward = []
    for t in t_list:
        k = list(tris[t]).index(u)
        forward.append(tris[t][(k + 1) % 3] == v)
    options = []
    for i in range(len(t_list)):
        for j in range(len(t_list)):
            if i != j and ok[i] and ok[j] and forward[i] != forward[j]:
                g = (turn[i] * (theta[j] - theta[i])) % (2 * np.pi)
                if 1e-4 < g < 2 * np.pi - 1e-4:
                    options.append((g, i, j))
    pairs, used = [], set()
    for g, i, j in sorted(options):          # tightest turns first, each triangle once
        if i not in used and j not in used:
            used |= {i, j}
            pairs.append((t_list[i], t_list[j]))
    return pairs


def split_bodies(pts, tris):
    """Groups of triangles to convert separately: one group when the mesh has no edge
    shared by more than two triangles (the usual case), otherwise bodies that touch go in
    different groups. Stray slivers (open or flat pieces of a few triangles) are dropped."""
    # (zero-area slivers are kept: along a T-junction they are what closes the gap)
    tris = tris[(tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 2] != tris[:, 0])]
    e = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(tris)), 3)
    key = e[:, 0].astype(np.int64) * (len(pts) + 1) + e[:, 1]
    order = np.argsort(key, kind="stable")
    _, start, cnt = np.unique(key[order], return_index=True, return_counts=True)
    if cnt.max(initial=0) <= 2:
        return [tris]

    parent, find = _union_find(len(tris))

    def join(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    touching = []
    for s, c in zip(start, cnt):
        ts = owner[order[s:s + c]]
        if c == 2:
            join(ts[0], ts[1])
        elif c > 2:
            u, v = e[order[s]]
            for a, b in _partners(pts, tris, list(ts), u, v):
                join(a, b)
            touching.append(ts)
    roots = np.array([find(i) for i in range(len(tris))])

    # keep closed pieces with real volume
    a, b, c = (pts[tris[:, k]] for k in range(3))
    tet = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6
    comps = {}
    for r in np.unique(roots):
        idx = np.nonzero(roots == r)[0]
        comps[r] = (idx, abs(tet[idx].sum()))
    biggest = max(v for _, v in comps.values())
    comps = {r: idx for r, (idx, v) in comps.items() if len(idx) >= 4 and v > 1e-6 * biggest}

    # bodies that meet at an edge must not be sewn together: colour them apart
    clash = {r: set() for r in comps}
    for ts in touching:
        rs = {roots[t] for t in ts} & comps.keys()
        for r in rs:
            clash[r] |= rs - {r}
    colour = {}
    for r in sorted(comps, key=lambda r: -len(comps[r])):
        used = {colour[q] for q in clash[r] if q in colour}
        colour[r] = next(k for k in range(len(comps)) if k not in used)
    groups = [np.concatenate([comps[r] for r in comps if colour[r] == k])
              for k in sorted(set(colour.values()))]
    return [tris[g] for g in groups]
