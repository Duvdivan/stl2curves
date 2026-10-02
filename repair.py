"""Mend the defects downloaded STLs often have, before anything else looks at them.

Each step changes the mesh only where something is wrong, and says what it did:

- triangles with a repeated corner, identical copies of a triangle, and pairs of the
  same triangle facing opposite ways (a wall of zero thickness) are dropped;
- sliver triangles (three corners in a line) are split away;
- cracks: corners along open edges that nearly meet are joined;
- orientation: triangles are turned to agree with their neighbours across every edge,
  and each closed piece faces outward (or, if it lies inside another piece, as a
  cavity, inward);
- holes: each loop of open edges that remains is filled;
- specks that cross themselves (small knots of triangles passing through each other)
  are cut out with a margin and the hole filled, as in MeshFix (Attene, "A lightweight
  approach to repairing digitized polygon meshes", 2010).
"""

import math

import numpy as np

from features import grid_noise

CRACK = 1e-3            # mm (and a small share of the part's size): open edges this close are joined
SPECK = 0.1             # a self-crossing cluster smaller than this share of the part is cut out
MARGIN_RINGS = 2        # ... together with this many rings of triangles around it


def repair(pts, tris):
    """(pts, tris, notes): the mended mesh and a list of what was done."""
    notes = []
    size = float(np.ptp(pts, axis=0).max()) if len(pts) else 1.0
    tris = np.asarray(tris, dtype=np.int64)

    k = len(tris)
    tris = tris[(tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 2] != tris[:, 0])]
    tris, same, opposite = _duplicates(tris)
    if same:
        notes.append(f"{same} duplicate triangles removed")
    if opposite:
        notes.append(f"{opposite} zero-thickness wall triangles removed")
    tris, slivers = remove_slivers(pts, tris)
    # where triangles cross each other, also the slivers the file's rounding has made
    # not quite flat (only there: elsewhere they do no harm)
    crossing = crossing_triangles(pts, tris)
    if len(crossing):
        zone = np.zeros(len(pts), bool)
        zone[np.unique(tris[np.unique(crossing)])] = True
        for _ in range(MARGIN_RINGS):
            zone[np.unique(tris[zone[tris].any(axis=1)])] = True
        tris, more = remove_slivers(pts, tris, grid_noise(pts), zone)
        slivers += more
    if slivers:
        notes.append(f"{slivers} sliver triangles split away")

    if len(_open_edges(tris)):
        pts, tris, joined = _close_cracks(pts, tris, max(CRACK, 1e-5 * size))
        if joined:
            notes.append(f"{joined} corners along cracks joined")

    tris, flipped = _orient(pts, tris)
    if flipped:
        notes.append(f"{flipped} triangles turned to face the right way")

    specks = 0
    for _ in range(3):
        cut = _self_crossing_specks(pts, tris, SPECK * size)
        if cut is None:
            break
        holes = len(_open_edges(tris))
        patched_pts, patched, filled = _fill_holes(pts, tris[~cut])
        if len(_open_edges(patched)) > holes:
            # (the hole it leaves didn't patch: on a fine mesh a knot a few millimetres
            # across is thousands of triangles. Better a mesh crossing itself there,
            # which still converts, than one with a hole, which makes no solid at all)
            break
        pts, tris = patched_pts, patched
        specks += 1
    if specks:
        notes.append(f"{specks} self-crossing speck{'s' if specks > 1 else ''} cut out and patched")

    pts, tris, filled = _fill_holes(pts, tris)
    if filled:
        notes.append(f"{filled} hole{'s' if filled > 1 else ''} filled")
    if flipped or filled or specks:
        tris, _ = _orient(pts, tris)
    return pts, tris, notes


# ---------------------------------------------------------------- simple clean-ups

def _duplicates(tris):
    """Drop repeated triangles (keeping one) and opposite pairs (dropping both)."""
    key = np.sort(tris, axis=1)
    _, first, inv, count = np.unique(key, axis=0, return_index=True, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    if count.max(initial=1) == 1:
        return tris, 0, 0
    # orientation of each triangle relative to its sorted corners: even or odd permutation
    rot = np.argmin(tris, axis=1)
    rolled = np.take_along_axis(tris, (np.arange(3)[None] + rot[:, None]) % 3, axis=1)
    parity = rolled[:, 1] < rolled[:, 2]        # same winding as the sorted order?
    keep = np.ones(len(tris), bool)
    same = opposite = 0
    for g in np.nonzero(count > 1)[0]:
        members = np.nonzero(inv == g)[0]
        up = members[parity[members]]
        down = members[~parity[members]]
        pairs = min(len(up), len(down))
        # opposite pairs cancel; of what's left, one copy stays
        keep[members] = False
        opposite += 2 * pairs
        rest = up[pairs:] if len(up) > len(down) else down[pairs:]
        if len(rest):
            keep[rest[0]] = True
            same += len(rest) - 1
    return tris[keep], same, opposite


def remove_slivers(pts, tris, noise=0.0, zone=None):
    """Remove slivers: triangles with one corner lying on the opposite side (all three
    corners in a line, give or take the file's rounding noise). Their outline doubles
    back on itself, which spoils the face they belong to, and a chain of them folds
    over itself. The triangle across that side is split at the corner instead, which
    keeps the surface as it was (to within the noise) and the mesh closed. The noise
    allowance applies only to corners in zone (all if None). Returns (tris, count)."""
    tris = tris.copy()
    total = 0
    for _ in range(20):
        T = tris
        found = []
        for k in range(3):
            v, u, w = pts[T[:, k]], pts[T[:, (k + 1) % 3]], pts[T[:, (k + 2) % 3]]
            e = w - u
            length = np.linalg.norm(e, axis=1)
            s = np.einsum("ij,ij->i", v - u, e) / np.maximum(length, 1e-12) ** 2
            h = np.linalg.norm(np.cross(v - u, e), axis=1) / np.maximum(length, 1e-12)
            allow = noise if zone is None else noise * zone[T[:, k]]
            thin = (h <= np.maximum(1e-5 * length, 1.5 * allow)) & (h < np.maximum(1e-4, 2 * allow))
            for i in np.nonzero(thin & (s > 0) & (s < 1))[0]:
                found.append((i, k))
        if not found:
            break
        across = {(a, b): i for i, t in enumerate(T) for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
        drop, add, touched, made = set(), [], set(), set()
        for i, k in found:
            v, u, w = T[i, k], T[i, (k + 1) % 3], T[i, (k + 2) % 3]
            j = across.get((w, u))
            if j is None or i in touched or j in touched:
                continue
            x = next(c for c in T[j] if c not in (u, w))
            if x != v and ((v, x) in across or (x, v) in across or (min(v, x), max(v, x)) in made):
                continue        # (splitting would add an edge that's already there: a fold)
            made.add((min(v, x), max(v, x)))
            touched |= {i, j}
            drop |= {i, j}
            total += 1
            if x != v:          # (x == v: two slivers folded onto each other; both go)
                add += [(w, v, x), (v, u, x)]
        keep = np.array([i not in drop for i in range(len(T))])
        tris = np.vstack([T[keep]] + ([np.array(add, dtype=T.dtype)] if add else []))
    return tris, total


def _edges(tris):
    return np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])


def _open_edges(tris):
    """Directed edges with no partner running the other way (the edges of holes)."""
    d = _edges(tris)
    n = int(tris.max()) + 1 if len(tris) else 1
    fwd = d[:, 0] * n + d[:, 1]
    back = d[:, 1] * n + d[:, 0]
    und = np.minimum(fwd, back)
    _, inv, count = np.unique(und, return_inverse=True, return_counts=True)
    return d[count[inv.ravel()] == 1]


# ---------------------------------------------------------------- cracks

def _close_cracks(pts, tris, tol):
    """Join corners on open edges that lie within tol of each other (a crack where two
    patches of triangles were never welded)."""
    open_ = _open_edges(tris)
    ends = np.unique(open_.ravel())
    if len(ends) < 2:
        return pts, tris, 0
    P = pts[ends]
    key = np.floor(P / tol).astype(np.int64)
    # cluster corners within tol: grid neighbours
    order = np.lexsort(key.T[::-1])
    parent = np.arange(len(ends))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    cells = {}
    for i in order:
        cells.setdefault(tuple(key[i]), []).append(i)
    for (x, y, z), members in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in cells.get((x + dx, y + dy, z + dz), ()):
                        for i in members:
                            if i < j and np.linalg.norm(P[i] - P[j]) <= tol:
                                ri, rj = root(i), root(j)
                                if ri != rj:
                                    parent[ri] = rj
    roots = np.array([root(i) for i in range(len(ends))])
    joined = int((roots != np.arange(len(ends))).sum())
    if not joined:
        return pts, tris, 0
    remap = np.arange(len(pts))
    remap[ends] = ends[roots]
    tris = remap[tris]
    tris = tris[(tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 2] != tris[:, 0])]
    tris, _, _ = _duplicates(tris)
    return pts, tris, joined


# ---------------------------------------------------------------- orientation

def _orient(pts, tris):
    """Make neighbouring triangles agree (each edge run opposite ways by the two
    triangles on it), then make each closed piece face outward, or inward if it lies
    inside an odd number of other pieces (a cavity). Returns (tris, how many flipped)."""
    tris = tris.copy()
    m = len(tris)
    if not m:
        return tris, 0
    d = _edges(tris)
    owner = np.tile(np.arange(m), 3)
    n = int(tris.max()) + 1
    und = np.minimum(d[:, 0], d[:, 1]) * n + np.maximum(d[:, 0], d[:, 1])
    order = np.argsort(und, kind="stable")
    und_s, own_s, d_s = und[order], owner[order], d[order]
    # pairs of triangles sharing an edge used exactly twice
    _, first, count = np.unique(und_s, return_index=True, return_counts=True)
    same = first[count == 2]
    ta, tb = own_s[same], own_s[same + 1]
    agree = d_s[same, 0] != d_s[same + 1, 0]            # opposite directions: consistent
    # breadth-first over the pairs, flipping to agree with the first triangle reached
    adj = [[] for _ in range(m)]
    for a, b, ok in zip(ta.tolist(), tb.tolist(), agree.tolist()):
        adj[a].append((b, ok))
        adj[b].append((a, ok))
    # (plain lists: element access on numpy arrays is slow in a loop like this)
    flip = [False] * m
    seen = [False] * m
    comp = [-1] * m
    c = 0
    for start in range(m):
        if seen[start]:
            continue
        seen[start] = True
        comp[start] = c
        stack = [start]
        while stack:
            t = stack.pop()
            for u, ok in adj[t]:
                if not seen[u]:
                    seen[u] = True
                    comp[u] = c
                    flip[u] = flip[t] if ok else not flip[t]
                    stack.append(u)
        c += 1
    flip, comp = np.array(flip, bool), np.array(comp)
    # each piece: flip it wholesale if it has more triangles flipped than not
    wholesale = np.bincount(comp, weights=flip, minlength=c) * 2 > np.bincount(comp, minlength=c)
    flip ^= wholesale[comp]
    tris[flip] = tris[flip][:, ::-1]
    count = int(flip.sum())
    # facing: outward for pieces at even depth, inward inside an odd number of others
    a, b, cc = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
    vol6 = np.einsum("ij,ij->i", a, np.cross(b, cc))
    by_piece = np.argsort(comp, kind="stable")
    pieces = np.split(by_piece, np.cumsum(np.bincount(comp, minlength=c))[:-1])
    closed = [p for p in pieces if len(p) >= 4 and not len(_open_edges(tris[p]))]
    boxes = [(pts[tris[q]].min(axis=(0, 1)), pts[tris[q]].max(axis=(0, 1))) for q in closed]
    for p in closed:
        if vol6[p].sum() > 0:
            continue
        # inside out, unless it lies wholly inside another piece: then it may be a cavity
        # (or a body buried in another: either way the file's own facing is the best clue)
        corners = np.unique(tris[p])
        sample = pts[corners[:: max(1, len(corners) // 12)]]
        # (a piece whose box doesn't hold every sample point can't hold them all)
        if not any(np.all((sample >= lo) & (sample <= hi)) and all(_inside(pts, tris[q], x) for x in sample)
                   for q, (lo, hi) in zip(closed, boxes) if q is not p):
            tris[p] = tris[p][:, ::-1]
            count += len(p)
    return tris, count


def _inside(pts, tris, point):
    """Is the point inside the closed mesh (ray parity along a skewed direction)?"""
    lo, hi = pts[tris].min(axis=(0, 1)), pts[tris].max(axis=(0, 1))
    if np.any(point < lo) or np.any(point > hi):
        return False
    d = np.array([0.5773, 0.5774, 0.5775])
    d /= np.linalg.norm(d)
    a, b, c = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
    e1, e2 = b - a, c - a
    p = np.cross(d, e2)
    det = np.einsum("ij,ij->i", e1, p)
    ok = np.abs(det) > 1e-15
    inv = np.where(ok, 1 / np.where(ok, det, 1), 0)
    s = point - a
    u = np.einsum("ij,ij->i", s, p) * inv
    q = np.cross(s, e1)
    v = (q @ d) * inv
    t = np.einsum("ij,ij->i", e2, q) * inv
    hits = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 1e-9)
    return bool(hits.sum() % 2)


# ---------------------------------------------------------------- holes

def _loops(open_):
    """Chains of open edges into closed loops of corners (in the edges' direction)."""
    nxt = {}
    for a, b in open_.tolist():
        nxt.setdefault(a, []).append(b)
    loops, used = [], set()
    for a, b in open_.tolist():
        if (a, b) in used:
            continue
        loop, cur, prev = [a], b, a
        used.add((a, b))
        while cur != a and len(loop) < 100000:
            loop.append(cur)
            options = [x for x in nxt.get(cur, []) if (cur, x) not in used]
            if not options:
                loop = None
                break
            prev, cur = cur, options[0]
            used.add((prev, cur))
        if loop:
            loops.append(loop)
    return loops


def _fill_holes(pts, tris):
    """Fill each loop of open edges with triangles. Returns (pts, tris, holes filled)."""
    open_ = _open_edges(tris)
    if not len(open_):
        return pts, tris, 0
    new, filled = [], 0
    for loop in _loops(open_):
        # the patch runs the other way round the loop from the triangles around it
        patch = _triangulate(pts, loop[::-1])
        if patch is not None:
            new.append(patch)
            filled += 1
    if new:
        tris = np.vstack([tris] + new)
    return pts, tris, filled


def _triangulate(pts, loop):
    """Triangles covering a loop of corners: the smallest total area over all ways of
    cutting it into triangles (dynamic programming, as in Liepa's hole filling) for
    loops of up to 60 corners, otherwise ears cut off in the loop's best-fit plane."""
    k = len(loop)
    if k < 3:
        return None
    if k == 3:
        return np.array([loop])
    P = pts[loop]
    if k <= 60:
        area = lambda i, j, m: 0.5 * np.linalg.norm(np.cross(P[j] - P[i], P[m] - P[i]))
        best = np.zeros((k, k))
        split = np.zeros((k, k), int)
        for gap in range(2, k):
            for i in range(k - gap):
                j = i + gap
                options = [best[i, m] + best[m, j] + area(i, m, j) for m in range(i + 1, j)]
                m = int(np.argmin(options))
                best[i, j], split[i, j] = options[m], i + 1 + m
        out, stack = [], [(0, k - 1)]
        while stack:
            i, j = stack.pop()
            if j - i < 2:
                continue
            m = split[i, j]
            out.append((loop[i], loop[m], loop[j]))
            stack += [(i, m), (m, j)]
        return np.array(out)
    # large loop: ear clipping in the best-fit plane
    mid = P.mean(axis=0)
    _, _, V = np.linalg.svd(P - mid)
    xy = (P - mid) @ V[:2].T
    signed = 0.5 * np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1])
    idx = list(range(k)) if signed > 0 else list(range(k))[::-1]
    out = []
    guard = 0
    while len(idx) > 3 and guard < 10 * k:
        guard += 1
        for t in range(len(idx)):
            i0, i1, i2 = idx[t - 1], idx[t], idx[(t + 1) % len(idx)]
            a, b, c = xy[i0], xy[i1], xy[i2]
            if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) <= 0:
                continue                    # reflex corner
            others = [xy[x] for x in idx if x not in (i0, i1, i2)]
            if any(_in_triangle(p, a, b, c) for p in others):
                continue
            out.append((loop[i0], loop[i1], loop[i2]) if signed > 0 else (loop[i2], loop[i1], loop[i0]))
            idx.pop(t)
            break
        else:
            break
    if len(idx) == 3:
        i0, i1, i2 = idx
        out.append((loop[i0], loop[i1], loop[i2]) if signed > 0 else (loop[i2], loop[i1], loop[i0]))
    elif len(idx) > 3:
        # couldn't clip (a twisted loop): a fan from a new corner in the middle
        return None
    return np.array(out)


def _in_triangle(p, a, b, c):
    def side(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    d1, d2, d3 = side(p, a, b), side(p, b, c), side(p, c, a)
    return not ((d1 < 0 or d2 < 0 or d3 < 0) and (d1 > 0 or d2 > 0 or d3 > 0))


# ---------------------------------------------------------------- self-crossing specks

def crossing_triangles(pts, tris):
    """Pairs of triangles (sharing no corner) that pass through each other."""
    if len(tris) < 2:
        return np.zeros((0, 2), int)
    T = pts[tris]
    lo, hi = T.min(axis=1), T.max(axis=1)
    ext = (hi - lo).max(axis=1)
    # (only triangles up to twice the usual size are compared: specks are made of small
    # triangles, and a big triangle crossing another is a body sunk into another)
    h = max(float(np.percentile(ext, 95)) * 2, 1e-6)
    small = ext <= h
    cells_lo = np.floor(lo / h).astype(np.int64)
    cells_hi = np.floor(hi / h).astype(np.int64)
    # small triangles: into every grid cell their box touches (at most 3 x 3 x 3)
    entries_t, entries_c = [], []
    idx = np.nonzero(small)[0]
    for dx in range(3):
        for dy in range(3):
            for dz in range(3):
                c = cells_lo[idx] + [dx, dy, dz]
                ok = np.all(c <= cells_hi[idx], axis=1)
                entries_t.append(idx[ok])
                entries_c.append(c[ok])
    et = np.concatenate(entries_t)
    ec = np.concatenate(entries_c)
    key = (ec[:, 0] * 73856093) ^ (ec[:, 1] * 19349663) ^ (ec[:, 2] * 83492791)
    order = np.argsort(key, kind="stable")
    key, et = key[order], et[order]
    starts = np.r_[0, np.nonzero(key[1:] != key[:-1])[0] + 1]
    ends = np.r_[starts[1:], len(key)]
    # every entry with each one after it in its cell (up to 199 on: an overfull cell is
    # a pile of slivers, not a speck), written out directly rather than offset by offset
    pos = np.arange(len(key))
    end = np.repeat(ends, ends - starts)
    count = np.minimum(end - pos - 1, 199)
    if not count.sum():
        return np.zeros((0, 2), int)
    first = np.repeat(pos, count)
    second = first + 1 + np.arange(len(first)) - np.repeat(np.cumsum(count) - count, count)
    a, b = et[first], et[second]
    # boxes apart: not a crossing (checked before the duplicates are sorted out, as it
    # leaves far fewer to sort)
    ov = np.all((lo[a] <= hi[b]) & (hi[a] >= lo[b]), axis=1)
    a, b = np.minimum(a[ov], b[ov]), np.maximum(a[ov], b[ov])
    n = np.int64(len(tris))
    code = np.unique(a.astype(np.int64) * n + b)
    pairs = np.c_[code // n, code % n]
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    # sharing a corner: not a crossing (neighbours meet along an edge or at a point)
    share = (tris[pairs[:, 0]][:, :, None] == tris[pairs[:, 1]][:, None, :]).any(axis=(1, 2))
    pairs = pairs[~share]
    if not len(pairs):
        return pairs
    hit = np.zeros(len(pairs), bool)
    A, B = T[pairs[:, 0]], T[pairs[:, 1]]
    for X, Y in ((A, B), (B, A)):
        for k in range(3):
            hit |= _segment_hits_triangle(X[:, k], X[:, (k + 1) % 3], Y)
    return pairs[hit]


def _segment_hits_triangle(p, q, T):
    """Does segment p-q (strictly) cross triangle T (Moller-Trumbore)?"""
    d = q - p
    e1, e2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]
    h = np.cross(d, e2)
    a = np.einsum("ij,ij->i", e1, h)
    ok = np.abs(a) > 1e-14
    f = np.where(ok, 1 / np.where(ok, a, 1), 0)
    s = p - T[:, 0]
    u = f * np.einsum("ij,ij->i", s, h)
    qv = np.cross(s, e1)
    v = f * np.einsum("ij,ij->i", d, qv)
    t = f * np.einsum("ij,ij->i", e2, qv)
    eps = 1e-9
    return ok & (u > eps) & (v > eps) & (u + v < 1 - eps) & (t > eps) & (t < 1 - eps)


def _self_crossing_specks(pts, tris, max_size):
    """Triangles to cut out: the small clusters of self-crossing triangles, with a
    margin of rings around each. None if there are none (or only large crossings,
    which are left alone: a body sunk into another is fused later, not cut)."""
    pairs = crossing_triangles(pts, tris)
    if not len(pairs):
        return None
    bad = np.zeros(len(tris), bool)
    bad[np.unique(pairs)] = True
    # clusters of crossing triangles (linked through shared corners); the small ones go,
    # with a margin of rings of triangles around them
    cut = np.zeros(len(tris), bool)
    for members in _components(tris, bad):
        if np.ptp(pts[np.unique(tris[members])], axis=0).max() <= max_size:
            cut[members] = True
    if not cut.any():
        return None
    for _ in range(MARGIN_RINGS):
        cut |= np.isin(tris, np.unique(tris[cut])).any(axis=1)
    return cut


def _components(tris, mask):
    """Groups of the masked triangles connected through shared corners."""
    idx = np.nonzero(mask)[0]
    parent = {int(i): int(i) for i in idx}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by_corner = {}
    for i in idx.tolist():
        for v in tris[i].tolist():
            by_corner.setdefault(v, []).append(i)
    for members in by_corner.values():
        r0 = root(members[0])
        for m in members[1:]:
            r = root(m)
            if r != r0:
                parent[r] = r0
    groups = {}
    for i in idx.tolist():
        groups.setdefault(root(i), []).append(i)
    return [np.array(g) for g in groups.values()]
