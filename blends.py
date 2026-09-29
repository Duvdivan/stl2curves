"""
Smooth blends: curved areas that no cylinder, cone, sphere or torus explains.

Where fillets meet at a corner, or a fillet follows a spline-shaped edge, the mesh
has a dense patch of small facets that detection either leaves flat or covers with
small mismatched pieces. Each such area becomes one smooth freeform face instead:
a surface fitted through the patch's mesh corners, kept tangent to the faces it
rolls into, and trimmed to the patch's outline (shared exactly with its neighbours).
"""
import math
from dataclasses import replace

import numpy as np

from features import Feature, Revolved, Sphere

SMOOTH_DEG = 50       # facets bending less than this meet smoothly
SPREAD_DEG = 40       # one blend face turns at most this far from its first facet and its mean direction
TANGENT_DEG = 35      # a neighbour this close in direction is rolled into tangentially
MIN_BEND_DEG = 2      # an unrecognised facet must turn at least this far against a neighbour
CORNER_SIZE = 1.5     # mm: a blend no bigger across than this is a corner...
CORNER_SPREAD_DEG = 85  # ...and may turn this far
MAX_SPLITS = 2        # a blend whose face won't fit is cut in two, at most this many times over
SMALL_FACET = 0.02    # unrecognised facets smaller than this share of the biggest can be blends
MAX_DEVIATION = 0.02  # a blend face must pass this close to every mesh corner (mm)
MAX_BULGE = 0.15      # ...and bow away from a facet by at most this share of its size
MAX_EDGE_GAP = 0.01   # ...and follow its outline edges this closely (mm)


class Blend:
    """Stand-in model for a blend: the surface itself is fitted when the face is built."""
    kind = "blend"
    line = circle = None

    def __init__(self, normal):
        self.normal_hint = normal


def _bend(mesh, a, b):
    return math.degrees(math.acos(float(np.clip(mesh.fn[a] @ mesh.fn[b], -1, 1))))


def _radius(model):
    if isinstance(model, Sphere):
        return model.r
    if isinstance(model, Revolved):
        if model.line:
            return model.line[0] if model.line[1] == 0 else None
        return model.circle[2]
    return None


def _scrap(mesh, f):
    """A small piece standing in for part of a surface it doesn't really fit: few facets,
    turning through a small angle, small for its radius."""
    r = _radius(f.model)
    if r is None or f.kind in ("ball", "wedge") or len(f.facets) >= 40:
        return False
    N = mesh.fn[f.facets]
    turn = math.degrees(math.acos(float(np.clip((N @ N.T).min(), -1, 1))))
    return turn < 45 and float(mesh.farea[f.facets].sum()) < r * r


def _ball_chains(mesh, features):
    """Sphere pieces of one radius touching one another but centred in different places:
    a rolling-ball fillet whose path isn't a circle (each short stretch fits a ball of the
    fillet radius, the balls' centres wander along the path). A real ball corner is one
    sphere, so three or more such pieces are all stand-ins."""
    balls = [k for k, f in enumerate(features) if isinstance(f.model, Sphere) and f.kind != "ball"]
    verts = {k: set(np.concatenate([mesh.fverts[x] for x in features[k].facets]).tolist()) for k in balls}
    link = {k: [] for k in balls}
    for i, a in enumerate(balls):
        for b in balls[i + 1:]:
            ma, mb = features[a].model, features[b].model
            if (abs(ma.r - mb.r) <= 0.05 * ma.r and np.linalg.norm(ma.c - mb.c) > 0.05 * ma.r
                    and verts[a] & verts[b]):
                link[a].append(b)
                link[b].append(a)
    out, seen = set(), set()
    for k in balls:
        if k in seen:
            continue
        group, stack = [], [k]
        seen.add(k)
        while stack:
            a = stack.pop()
            group.append(a)
            for b in link[a]:
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
        if len(group) >= 3:
            out.update(group)
    return out


def add_blends(mesh, features):
    """Features with small mismatched pieces and unrecognised curved facets regrouped into
    smooth blends. Each blend remembers the pieces it replaced (`parts`), so they can be
    put back if its face can't be built."""
    nf = len(mesh.farea)
    owner = np.full(nf, -1)
    for k, f in enumerate(features):
        owner[f.facets] = k
    scraps = {k for k, f in enumerate(features) if _scrap(mesh, f)} | _ball_chains(mesh, features)
    small = SMALL_FACET * mesh.farea.max()

    def candidate(a):
        if owner[a] >= 0:
            return owner[a] in scraps
        # a facet that barely turns against any neighbour belongs to a flat face, not a blend
        bends = [_bend(mesh, a, b) for b in mesh.nbrs[a]]
        return mesh.farea[a] <= small and any(MIN_BEND_DEG <= x < SMOOTH_DEG for x in bends)

    cand = np.array([candidate(a) for a in range(nf)])

    # Regions grow a unit at a time: an unrecognised facet, or a whole piece (so a piece
    # can be put back whole if its blend fails). Every facet must stay within SPREAD_DEG
    # of the region's first facet and of its mean direction: no further than one smooth
    # face can follow. A region with no curve in it (flat facets only) stays as it is.
    def unit(a):
        return [int(x) for x in features[owner[a]].facets] if owner[a] >= 0 else [int(a)]

    def angle(n, m):
        return math.degrees(math.acos(float(np.clip(n @ m, -1, 1))))

    blends = []
    free = set(np.nonzero(cand)[0].tolist())
    while free:
        seed_unit = unit(max(free, key=lambda a: mesh.farea[a]))
        seed_n = mesh.fn[seed_unit].T @ mesh.farea[seed_unit]
        seed_n = seed_n / np.linalg.norm(seed_n)
        region, stack, normal = list(seed_unit), list(seed_unit), seed_n * mesh.farea[seed_unit].sum()
        free.difference_update(seed_unit)
        while stack:
            a = stack.pop()
            for b in mesh.nbrs[a]:
                if b not in free or _bend(mesh, a, b) >= SMOOTH_DEG:
                    continue
                u = unit(b)
                if any(x not in free for x in u):
                    continue
                mean = normal / np.linalg.norm(normal)
                turn = max(max(angle(mesh.fn[x], seed_n), angle(mesh.fn[x], mean)) for x in u)
                if turn > SPREAD_DEG:
                    # a small corner (where rounded edges meet) may turn further: one
                    # N-sided patch spans it, bounded by the exact faces around it
                    if turn > CORNER_SPREAD_DEG:
                        continue
                    P = mesh.pts[np.concatenate([mesh.fverts[x] for x in region + u])]
                    if np.linalg.norm(np.ptp(P, axis=0)) > CORNER_SIZE:
                        continue
                free.difference_update(u)
                region += u
                stack += u
                normal = normal + mesh.fn[u].T @ mesh.farea[u]
        members = set(region)
        inner_bends = [_bend(mesh, a, b) for a in region for b in mesh.nbrs[a] if b in members]
        if len(region) < 3 or max(inner_bends, default=0) < 0.5:
            continue
        parts = sorted({int(owner[a]) for a in region if owner[a] >= 0})
        blends.append((np.array(sorted(region)), parts, normal / np.linalg.norm(normal)))

    out, dropped = [], set()
    for facets, parts, normal in blends:
        N = mesh.fn[facets]
        turn = math.degrees(math.acos(float(np.clip((N @ N.T).min(), -1, 1))))
        area = float(mesh.farea[facets].sum())
        out.append(_blend(mesh, facets, tuple(features[k] for k in parts)))
        dropped.update(parts)
    return [f for k, f in enumerate(features) if k not in dropped] + out


def _blend(mesh, facets, parts, depth=0):
    N = mesh.fn[facets]
    turn = math.degrees(math.acos(float(np.clip((N @ N.T).min(), -1, 1))))
    normal = N.T @ mesh.farea[facets]
    area = float(mesh.farea[facets].sum())
    return Feature(model=Blend(normal / np.linalg.norm(normal)), label="smooth blend",
                   detail=f"{len(facets)} facets, turns {turn:.0f} deg", convex=True, kind="blend",
                   facets=np.asarray(facets), change=0.0, tolerance=0.02 * area, worst=0.0,
                   parts=parts, depth=depth)


def split(mesh, blend):
    """A blend cut in two across its longest direction (each piece it replaced going
    whole to the half holding most of it), or None if it can't be cut any further."""
    if blend.depth >= MAX_SPLITS or len(blend.facets) < 8:
        return None
    C = mesh.fcent[blend.facets]
    axis = np.linalg.svd(C - C.mean(axis=0))[2][0]
    t = (C - C.mean(axis=0)) @ axis
    cut = np.median(t)
    side = {int(f): bool(x > cut) for f, x in zip(blend.facets, t)}
    halves = [[], []], [[], []]          # (facets, parts) for each side
    in_part = set()
    for p in blend.parts:
        votes = [side[int(f)] for f in p.facets if int(f) in side]
        h = int(sum(votes) * 2 > len(votes))
        halves[h][0].extend(int(f) for f in p.facets)
        halves[h][1].append(p)
        in_part.update(int(f) for f in p.facets)
    for f in blend.facets:
        if int(f) not in in_part:
            halves[int(side[int(f)])][0].append(int(f))
    out = []
    for facets, parts in halves:
        if len(facets) < 3:
            return None
        out.append(_blend(mesh, np.array(sorted(facets)), tuple(parts), blend.depth + 1))
    return out
