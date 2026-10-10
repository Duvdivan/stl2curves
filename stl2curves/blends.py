"""
Smooth blends: curved areas that no cylinder, cone, sphere or torus explains.

Where fillets meet at a corner, or a fillet follows a spline-shaped edge, the mesh
has a dense patch of small facets that detection either leaves flat or covers with
small mismatched pieces. Each such area becomes one smooth freeform face instead:
an N-sided patch spanning the area's outline (shared exactly with its neighbours) and
passing through its mesh corners.
"""
import math
from collections import deque

import numpy as np

from . import extrude
from . import freeform
from . import pipes
from .accuracy import SCALE as ACCURACY
from .features import Feature, Revolved, Sphere

SMOOTH_DEG = 50       # facets bending less than this meet smoothly
SPREAD_DEG = 40       # one blend face turns at most this far from its first facet and its mean direction
TANGENT_DEG = 35      # a neighbour this close in direction is rolled into tangentially
MIN_BEND_DEG = 2      # an unrecognised facet must turn at least this far against a neighbour
MAX_FACETS = 60       # a bigger smooth area is covered by several blends (each a manageable fit)
CORNER_SIZE = 1.5     # mm: a blend no bigger across than this is a corner...
CORNER_SPREAD_DEG = 85  # ...and may turn this far
MAX_SPLITS = 2        # a blend whose face won't fit is cut in two, at most this many times over
SMALL_FACET = 0.02    # unrecognised facets smaller than this share of the biggest can be blends
MAX_DEVIATION = 0.02 * ACCURACY  # a blend face must pass this close to every mesh corner (mm)
MAX_BULGE = 0.15      # ...and bow away from a facet by at most this share of its size
MAX_EDGE_GAP = 0.01 * ACCURACY   # ...and follow its outline edges this closely (mm)
CREASE_DEG = 20       # a freeform area (freeform.py) doesn't reach across a sharper bend than this
FREEFORM_TURN_DEG = 3  # ...and turns at least this far from its mean direction somewhere
FREEFORM_MIN_RADIUS = 2.0   # mm: a piece rounder than this (a small rounded edge) stays out of one
FREEFORM_SPLITS = 4   # an area no surface fits is halved by facet direction at most this many times over
STRIP_DEG = 45        # an exact patch turning less than this beside freeform pieces may join them
FREEFORM_MIN_SIZE = 5.0     # mm across: a smaller area of MAX_FACETS or fewer is left to the blends
GROW_COS = math.cos(math.radians(10))   # a facet a freeform piece takes in faces within 10 deg of it
CARVE_MIN = 200       # facets: a freeform patch this big blamed for trouble loses just its facets near it


class Blend:
    """Stand-in model for a blend: the surface itself is fitted when the face is built."""
    kind = "blend"
    line = circle = None

    def __init__(self, normal):
        self.normal_hint = normal


def _bend(mesh, a, b):
    return math.degrees(math.acos(float(np.clip(mesh.fn[a] @ mesh.fn[b], -1, 1))))


def _radius(model, P=None):
    """The surface's radius (a cone's: its mean radius at the points P)."""
    if isinstance(model, Sphere):
        return model.r
    if isinstance(model, Revolved):
        if model.line:
            if model.line[1] == 0:
                return model.line[0]
            return float(model.local(P)[0].mean()) if P is not None else None
        return model.circle[2]
    return None


def _scrap(mesh, f, any_count=False):
    """A small piece standing in for part of a surface it doesn't really fit: few facets
    (or any number, any_count: on a finely meshed part a narrow strip has dozens),
    turning through a small angle, small for its radius."""
    r = _radius(f.model, mesh.fcent[f.facets])
    if r is None or f.kind in ("ball", "wedge") or (len(f.facets) >= 40 and not any_count):
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
    # (first the pieces of fillets along curved edges joined into pipes: they are pieces
    # of tori, cylinders and spheres the passes below would take for scraps)
    features = pipes.add_pipes(mesh, features)
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

    # a unit: an unrecognised facet, or a whole piece (so a piece can be put back whole)
    def unit(a):
        return [int(x) for x in features[owner[a]].facets] if owner[a] >= 0 else [int(a)]

    free = set(np.nonzero(cand)[0].tolist())
    out, dropped = [], set()
    # Big smooth areas first, each as one freeform surface if one fits; the rest (and
    # any that don't fit) as smooth blends of up to MAX_FACETS facets. (A freeform area
    # also takes small facets that bend only slightly against their neighbours: on a
    # finely meshed surface every bend is slight.)
    stand_ins = scraps | {k for k, f in enumerate(features) if _scrap(mesh, f, any_count=True)}
    stand_ins = {k for k in stand_ins
                 if (_radius(features[k].model, mesh.fcent[features[k].facets]) or 0) >= FREEFORM_MIN_RADIUS}
    soft = {a for a in range(nf) if (owner[a] in stand_ins if owner[a] >= 0
                                     else mesh.farea[a] <= small and mesh.nbrs[a])}
    # (profiles pushed along a design direction first: a freeform patch would take
    # their strips too, but as a surface that isn't one)
    ext, ext_parts, ext_taken = extrude.add_extrusions(mesh, features, soft, owner, unit)
    out += ext
    dropped.update(ext_parts)
    soft -= ext_taken
    free -= ext_taken
    pieces =[p for area in _areas(mesh, soft, unit) for p in _freeform_pieces(mesh, area, unit)]
    pieces = _absorb_strips(mesh, features, owner, pieces, ext_parts)
    pieces = _grown(mesh, pieces, soft, unit)
    for facets, model, dev in pieces:
        parts = sorted({int(owner[a]) for a in facets if owner[a] >= 0})
        out.append(freeform.feature(mesh, model, facets, tuple(features[k] for k in parts), dev))
        dropped.update(parts)
        free.difference_update(facets.tolist())
    for facets, parts in _regions(mesh, free, unit, lambda a: owner[a]):
        out.append(_blend(mesh, facets, tuple(features[k] for k in parts)))
        dropped.update(parts)
    return [f for k, f in enumerate(features) if k not in dropped] + out


def merge(mesh, features, facets):
    """(features, merged): an area someone painted (stl2curves.studio: one the conversion
    left as triangles, a fillet it missed) made one face. A feature they cover half or
    more of goes into it whole, one they barely touch stays as it was, and flat facets
    the brush caught (none of whose neighbours meet it at a gentle bend) stay flat. The
    area is one cylinder, cone or sphere if it all lies on one, else one freeform surface
    if one fits, else one smooth blend (cut in two by the build if no fill fits); either
    way it gives back the features it took if its face can't be built. merged is None if
    there is nothing to merge."""
    from .features import Region, _whole_region
    chosen = {int(a) for a in facets if len(mesh.nbrs[int(a)])}
    keep, taken = [], []
    for f in features:
        inside = sum(int(a) in chosen for a in f.facets)
        if inside and 2 * inside >= len(f.facets):
            taken.append(f)
        else:
            keep.append(f)
            chosen.difference_update(int(a) for a in f.facets)
    for f in taken:
        chosen.update(int(a) for a in f.facets)
    if len(chosen) < 2:
        return features, None
    area = np.array(sorted(chosen))
    merged = None
    try:
        best = _whole_region(mesh, Region(mesh, area))
    except Exception:
        best = None
    if best is not None and len(best[1]) == len(area):
        merged = best[0]
        merged.parts = tuple(taken)
    if merged is None:
        model, dev = freeform.fit(mesh, area)
        merged = (freeform.feature(mesh, model, area, tuple(taken), dev) if model is not None
                  else _blend(mesh, area, tuple(taken)))
    merged.detail = "painted: " + merged.detail
    return keep + [merged], merged


def fallback(mesh, feature):
    """Smooth blends of up to MAX_FACETS facets over a freeform patch's area (for when its
    face can't be built), the pieces it replaced kept whole."""
    owner = np.full(len(mesh.farea), -1)
    for k, p in enumerate(feature.parts):
        owner[p.facets] = k

    def unit(a):
        return [int(x) for x in feature.parts[owner[a]].facets] if owner[a] >= 0 else [int(a)]

    free = set(int(a) for a in feature.facets)
    regions = _regions(mesh, free, unit, lambda a: owner[a])
    used = {k for _, parts in regions for k in parts}
    # (pieces no blend took, and the facets left over, stay as they were)
    return ([_blend(mesh, facets, tuple(feature.parts[k] for k in parts)) for facets, parts in regions]
            + [p for k, p in enumerate(feature.parts) if k not in used])


def pipe_fallback(mesh, feature):
    """The pieces a pipe joined, as they were, and smooth blends of up to MAX_FACETS
    facets over the loose facets it took in besides (for when its face can't be built)."""
    taken = {int(a) for p in feature.parts for a in p.facets}
    free = {int(a) for a in feature.facets} - taken
    regions = _regions(mesh, free, lambda a: [a], lambda a: -1)
    return list(feature.parts) + [_blend(mesh, facets, ()) for facets, _ in regions]


def _freeform_pieces(mesh, facets, unit, depth=0):
    """[(facets, surface, deviation)]: the area as one freeform patch if a surface fits
    it, else its two halves by facet direction (a U-shaped area turns too far for one
    height field; a crease runs between faces facing different ways), each in turn.
    Pieces too small for a freeform patch are left to the blends."""
    if len(facets) < 3 or _turn(mesh, facets) < FREEFORM_TURN_DEG:
        return []
    if len(facets) <= MAX_FACETS:
        # (a few long facets can span a big area too)
        P = mesh.pts[np.unique(np.concatenate([mesh.fverts[a] for a in facets]))]
        if np.linalg.norm(np.ptp(P, axis=0)) < FREEFORM_MIN_SIZE:
            return []
    model, dev = freeform.fit(mesh, facets)
    if model is not None:
        return [(facets, model, dev)]
    if depth >= FREEFORM_SPLITS:
        return []
    # (split square to the direction the facets' normals vary most, through their mean)
    N, w = mesh.fn[facets], mesh.farea[facets]
    mean = N.T @ w / w.sum()
    X = N - mean
    axis = np.linalg.eigh((X * w[:, None]).T @ X)[1][:, -1]
    side = dict(zip(facets.tolist(), (X @ axis > 0).tolist()))
    halves = {True: set(), False: set()}
    for a in facets.tolist():
        if a not in halves[True] and a not in halves[False]:
            u = unit(a)         # (a piece goes whole to the side most of it is on)
            halves[sum(side[x] for x in u) * 2 > len(u)].update(u)
    out = []
    for members in halves.values():
        for piece in _areas(mesh, members, unit):
            out += _freeform_pieces(mesh, piece, unit, depth + 1)
    return out


def _absorb_strips(mesh, features, owner, pieces, skip=()):
    """Freeform pieces with the narrow exact patches beside them taken in, where one
    surface fits both: a curvature-continuous fillet (a spline in the original) fits a
    cylinder exactly along a narrow band in its middle, with freeform pieces either side;
    it is one face. (Only patches turning less than STRIP_DEG, met smoothly.)"""
    if not pieces:
        return pieces
    piece_of = np.full(len(mesh.farea), -1)
    for i, (facets, _, _) in enumerate(pieces):
        piece_of[facets] = i
    pieces = list(pieces)
    limit = math.cos(math.radians(CREASE_DEG))
    for k, f in enumerate(features):
        m = f.model
        if k in skip:
            continue
        if f.kind not in ("revolve", "trimmed") or not isinstance(m, (Revolved, Sphere))                 or _turn(mesh, f.facets) * 2 >= STRIP_DEG or (piece_of[f.facets] >= 0).any():
            continue
        # the pieces it meets smoothly
        touching = {int(piece_of[b]) for a in f.facets for b in mesh.nbrs[a]
                    if piece_of[b] >= 0 and pieces[piece_of[b]] is not None and mesh.fn[a] @ mesh.fn[b] >= limit}
        if not touching:
            continue
        union = np.unique(np.concatenate([f.facets] + [pieces[i][0] for i in touching]))
        model, dev = freeform.fit(mesh, union)
        if model is None:
            continue
        for i in touching:
            pieces[i] = None
        pieces.append((union, model, dev))
        piece_of[union] = len(pieces) - 1
    return [p for p in pieces if p is not None]


def _grown(mesh, pieces, soft, unit):
    """Freeform pieces with the loose facets beside them taken in where those already lie
    on the piece's surface: an area halved by direction leaves fringes too small to fit
    on their own, and a row of triangles a hundredth of a millimetre wide left between
    two faces has edges shorter than any sewing tolerance."""
    taken = set()
    for facets, _, _ in pieces:
        taken.update(facets.tolist())
    out = []
    for facets, model, dev in pieces:
        members = set(facets.tolist())
        stack = list(members)
        while stack:
            a = stack.pop()
            for b in mesh.nbrs[a]:
                if b in members or b in taken or b not in soft:
                    continue
                u = unit(b)
                if any(x in taken or x not in soft for x in u):
                    continue
                P = mesh.pts[np.unique(np.concatenate([mesh.fverts[x] for x in u]))]
                if np.abs(model.signed(P)).max() > freeform.FIT_DEV or                         (np.einsum("ij,ij->i", mesh.fn[u], model.normal(mesh.fcent[u])) ** 2).min() < GROW_COS ** 2:
                    continue
                members.update(u)
                taken.update(u)
                stack += u
        out.append((np.array(sorted(members)), model, dev))
    return out


def _turn(mesh, facets):
    """How far (degrees) the facets turn from their mean direction at most."""
    n = mesh.fn[facets].T @ mesh.farea[facets]
    if np.linalg.norm(n) < 1e-9 * mesh.farea[facets].sum():
        return 180.0            # (facing every way: a whole ring, say)
    n = n / np.linalg.norm(n)
    return math.degrees(math.acos(float(np.clip((mesh.fn[facets] @ n).min(), -1, 1))))


def carve(mesh, feature, points, reach):
    """A big freeform patch blamed for trouble at these points (a gap at its outline,
    say) less its facets within reach of them, which are left flat, and less any piece
    it took in that reaches them, given back whole: [the patch, the pieces...]; None if
    the patch is small or would lose half of itself (then its blends are better)."""
    if points is None or not len(points) or len(feature.facets) < CARVE_MIN:
        return None
    from scipy.spatial import cKDTree
    facets = feature.facets
    size = np.array([np.linalg.norm(mesh.pts[mesh.fverts[a]] - mesh.fcent[a], axis=1).max() for a in facets])
    d, _ = cKDTree(np.asarray(points, float)).query(mesh.fcent[facets])
    cut = set(facets[d <= reach + size].tolist())
    back = [p for p in feature.parts if cut & set(p.facets.tolist())]
    for p in back:
        cut.update(p.facets.tolist())
    keep = np.array([a for a in facets.tolist() if a not in cut])
    if not cut or len(keep) < len(facets) / 2:
        return None
    parts = tuple(p for p in feature.parts if all(p is not b for b in back))
    return [freeform.feature(mesh, feature.model, keep, parts, feature.model.dev)] + back


def _areas(mesh, free, unit):
    """Smooth areas of the free facets (whole units), not reaching across a crease."""
    out, seen = [], set()
    limit = math.cos(math.radians(CREASE_DEG))
    for seed in sorted(free):
        if seed in seen:
            continue
        area = unit(seed)
        if any(x not in free for x in area):
            continue
        seen.update(area)
        stack = list(area)
        while stack:
            a = stack.pop()
            for b in mesh.nbrs[a]:
                if b in seen or b not in free or mesh.fn[a] @ mesh.fn[b] < limit:
                    continue
                u = unit(b)
                if any(x not in free or x in seen for x in u):
                    continue
                seen.update(u)
                area += u
                stack += u
        out.append(np.array(sorted(area)))
    return out


def _regions(mesh, free, unit, owner_of):
    """Blend regions over the free facets: [(facets, indices of the pieces in it)]."""
    # Regions grow a unit at a time. Every facet must stay within SPREAD_DEG of the
    # region's first facet and of its mean direction: no further than one smooth face
    # can follow. A region with no curve in it (flat facets only) stays as it is.
    def angle(n, m):
        return math.degrees(math.acos(float(np.clip(n @ m, -1, 1))))

    blends = []
    free = set(free)
    while free:
        seed_unit = unit(max(free, key=lambda a: mesh.farea[a]))
        seed_n = mesh.fn[seed_unit].T @ mesh.farea[seed_unit]
        seed_n = seed_n / np.linalg.norm(seed_n)
        region, normal = list(seed_unit), seed_n * mesh.farea[seed_unit].sum()
        stack = deque(seed_unit)        # breadth first: compact regions
        free.difference_update(seed_unit)
        while stack:
            a = stack.popleft()
            for b in mesh.nbrs[a]:
                if b not in free or _bend(mesh, a, b) >= SMOOTH_DEG:
                    continue
                u = unit(b)
                if any(x not in free for x in u) or len(region) + len(u) > MAX_FACETS:
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
        parts = sorted({int(owner_of(a)) for a in region if owner_of(a) >= 0})
        blends.append((np.array(sorted(region)), parts))
    return blends


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
