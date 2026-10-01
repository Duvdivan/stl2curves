"""
Fillets between flat faces, found the way a designer sees them.

Two flat faces that would meet at a sharp edge, rounded off with a fillet: the fillet
is a cylinder touching both faces (a ball of the fillet radius rolled along the
edge). Its axis runs parallel to the edge where the faces would meet, so given the two
faces the only unknown is the radius, and every mesh corner on the rounded band
between them pins it down. (Constant-radius "rolling ball" blend recovery between
primary surfaces, as in reverse-engineering work by Várady, Benkő, Kós and Martin.)

Fillets are everywhere in printed parts, and a mesh often keeps too little accuracy to
show them cleanly, so this is forgiving where the general surface passes can't be:

- The whole band between the two faces is collected first (facets turning steadily
  from one face's normal to the other's, lying between the faces), and the radius is
  the one most of its corners agree on, so a few stray corners don't spoil it.
- A band whose corners stay within FILLET_TOL of that fillet is taken: looser than the
  mesh's own tolerance, since the faces either side pin the surface down. The faces
  must run along the whole band (a narrow strip of a curved wall is no face for a
  rounded edge that follows the curve).
- One fillet radius is usually used all over a design. Bands that miss on their own are
  tried again with each radius found with confidence elsewhere in the part held fixed
  (a 4.99 or 4.9 mm look-alike next to clean 5 mm fillets is a 5 mm fillet), within
  PRIOR_TOL.

It runs before the passes that propose surfaces from pairs of neighbouring facets,
which would otherwise cut a coarsely meshed fillet into narrow strips of their own.
Finely cut pins aren't mistaken for fillets: the faces a fillet joins must be clearly
bigger than the facets of its band.
"""
import math

import numpy as np

import features as F
from features import Revolved

MIN_FACE_AREA = 2       # mm^2: flat facets at least this big can be the faces a fillet joins
MIN_ANGLE_DEG = 30      # the faces must meet at least this far from flat...
MAX_ANGLE_DEG = 170     # ...and not be (nearly) parallel
LINK_RINGS = 24         # a fillet's band is at most this many facets across
FACE_TO_STRIP = 3       # each face a fillet joins is at least this many times its biggest band facet
FILLET_TOL = 0.02       # mm: how far most band corners may lie off the fillet found from them ...
GROW_TOL = 0.05         # mm: ... and its facets' corners at most (meshes lose accuracy at fillets)
REL_TOL = 0.04          # ... but neither more than this share of its radius
SPAN_SLACK_DEG = 8      # deg: a fillet may turn this much more than its faces do
CONFIDENT_TOL = 0.005   # mm: fillets fitting this well (or cylinders and tori found strictly)
                        # vouch for their radius elsewhere in the part
PRIOR_RANGE = 0.05      # a band whose own radius is within this share of a vouched one ...
PRIOR_TOL = 0.05        # mm: ... is that radius if its corners lie this close to it
FACE_LENGTH = 1.0       # a fillet, and the faces it joins alongside each other, run this many radii along it
FACE_WIDTH = 0.3        # each face it joins is this many radii wide across it
ROUGH_AGREE = 0.7       # (faces stopping short of their meeting line by distances agreeing within e^this)
GROW_FACTOR = 3         # its facets may lie this many times as far off as its middle does
ON_SHARE = 75           # percent of a band's corners that must lie within the tolerance
CORNER_AGREE = 0.2      # a band facet's corners give radii agreeing within e^this ...
WEDGE_SLACK = 0.01      # mm: how far outside the faces' wedge a band corner may lie
DEPTH_AGREE = 0.3       # ... and agreeing with the rest of the band's within e^this


def _adjacency(mesh):
    """Every facet's neighbours across any shared edge, however sharp the bend."""
    nbrs = [set() for _ in range(len(mesh.farea))]
    for fs in mesh.edge_facets.values():
        for a in fs:
            for b in fs:
                if a != b:
                    nbrs[a].add(b)
    return nbrs


def find(mesh, regions, taken=None, known=()):
    """Fillet features between pairs of flat faces. Marks their facets as used in the
    regions so later passes leave them alone. known: radii found with confidence so far
    (cylinders, tori), which fillets can take on."""
    nbrs = _adjacency(mesh)
    used = np.zeros(len(mesh.farea), bool) if taken is None else taken.copy()
    faces = [f for f in np.argsort(-mesh.farea) if mesh.farea[f] >= MIN_FACE_AREA and not used[f]]
    plane_pt = {f: mesh.pts[mesh.fverts[f][0]] for f in faces}
    # Pairs of flat faces linked by a band of smaller facets (at most LINK_RINGS of
    # them across), meeting at a real angle; the biggest pairs first (real faces before
    # the strips of a fillet).
    is_face = np.zeros(len(mesh.farea), bool)
    is_face[faces] = True
    lo_c, hi_c = math.cos(math.radians(MAX_ANGLE_DEG)), math.cos(math.radians(MIN_ANGLE_DEG))
    pairs = set()
    for A in faces:
        small = mesh.farea[A] / FACE_TO_STRIP
        seen, frontier = {A}, [A]
        for _ in range(LINK_RINGS):
            step = []
            for g in frontier:
                for h in nbrs[g]:
                    if h in seen:
                        continue
                    seen.add(h)
                    if is_face[h] and lo_c <= mesh.fn[A] @ mesh.fn[h] <= hi_c:
                        pairs.add((min(A, h), max(A, h)))
                    if not used[h] and mesh.farea[h] <= small:
                        step.append(h)      # (a band facet: the far face is further on)
            frontier = step
            if not frontier:
                break
    pairs = sorted(((min(mesh.farea[A], mesh.farea[B]), A, B) for A, B in pairs), key=lambda x: -x[0])

    noise = getattr(mesh, "noise", 0.0)
    corners = _Corners(mesh)
    cache = np.full(len(mesh.farea), -1.0)       # (verdicts on facets, reset after each try)
    vouched = sorted(set(round(float(r), 6) for r in known if r > 0))
    found, waiting = [], []
    # First only the fillets that fit closely: they vouch for their radius. Then the
    # rest, with every radius vouched for in hand (a 2.52 mm band near clean 2.5 mm
    # fillets is a 2.5 mm fillet).
    for round_ in (0, 1):
        todo = pairs if round_ == 0 else waiting
        for item in todo:
            _, A, B = item
            if used[A] or used[B]:
                continue
            if round_ == 0 and not _could_join(mesh, A, B, plane_pt):
                continue
            got = None
            for X, Y in ((A, B), (B, A)):    # (from either face: one may border others too)
                got = _fillet(mesh, nbrs, used, X, Y, plane_pt, noise, vouched, corners, cache)
                if got is not None:
                    break
            if got is None or (round_ == 0 and got[5] > CONFIDENT_TOL + noise):
                if round_ == 0:
                    waiting.append(item)     # (taken in the second round, if at all)
                continue
            facets, sign, a, d, r, err = got
            feature = _feature(mesh, nbrs, A, B, a, d, r, sign, facets)
            if feature is None:
                continue
            used[facets] = True
            found.append(feature)
            if err <= CONFIDENT_TOL + noise:
                vouched = sorted(set(vouched) | {round(r, 6)})
    return found


def _could_join(mesh, A, B, plane_pt):
    """Could faces A and B be joined by a fillet, judged by the faces alone? A fillet of
    radius r touches each face r tan(theta/2) short of the line where they would meet
    (theta the angle between their normals): both faces stop about equally short of it,
    and they must be wide enough, and run alongside each other far enough, for that r.
    (A quick look, before the band between them is gathered: on a finely cut bore, every
    pair of its strips would otherwise be tried.)"""
    n1, n2 = mesh.fn[A], mesh.fn[B]
    d = np.cross(n1, n2)
    d /= np.linalg.norm(d)
    p = np.linalg.solve(np.array([n1, n2, d]), [n1 @ plane_pt[A], n2 @ plane_pt[B], 0.0])
    near, width, along = [], [], []
    for X, n in ((A, n1), (B, n2)):
        P = mesh.pts[mesh.fverts[X]] - p
        across = P @ np.cross(n, d)
        near.append(float(np.abs(across).min()) if (across > 0).all() or (across < 0).all() else 0.0)
        width.append(float(np.ptp(across)))
        along.append((float((P @ d).min()), float((P @ d).max())))
    theta = math.acos(max(-1.0, min(1.0, float(n1 @ n2))))
    if min(near) <= 0:
        return True                 # (a face reaching the line: a sharp edge, or nothing to tell)
    if abs(math.log(near[0] / near[1])) > ROUGH_AGREE:
        return False
    r = math.sqrt(near[0] * near[1]) / math.tan(theta / 2)
    run = min(along[0][1], along[1][1]) - max(along[0][0], along[1][0])
    rough = math.exp(-ROUGH_AGREE) * r
    return min(width) >= FACE_WIDTH * rough and run >= FACE_LENGTH * rough


class _Corners:
    """Every facet's corners in one array, to look at many facets at once."""

    def __init__(self, mesh):
        self.count = np.array([len(v) for v in mesh.fverts])
        self.start = np.r_[0, np.cumsum(self.count)[:-1]]
        self.ids = np.concatenate(mesh.fverts)
        self.pts = mesh.pts

    def of(self, gs):
        """(corner points of facets gs, where each facet's run of them starts)."""
        n = self.count[gs]
        first = np.r_[0, np.cumsum(n)[:-1]]
        k = np.repeat(self.start[gs] - first, n) + np.arange(int(n.sum()))
        return self.pts[self.ids[k]], first


def _fillet(mesh, nbrs, used, A, B, plane_pt, noise, vouched, corners, cache):
    """The rolling-ball fillet joining faces A and B, if there is one: (facets, side, a
    point on its axis, the axis direction, radius, worst corner gap), or None.

    The band between the faces is gathered loosely (facets whose corners each lie on
    some fillet between the faces, of about the same radius), its radius read from the
    corners where both faces run alongside (not where it meets a rounded corner at its
    end), and the fillet is then the band's facets that lie on that cylinder."""
    n1, n2 = mesh.fn[A], mesh.fn[B]
    d = np.cross(n1, n2)
    d /= np.linalg.norm(d)
    p = np.linalg.solve(np.array([n1, n2, d]), [n1 @ plane_pt[A], n2 @ plane_pt[B], 0.0])
    small = min(mesh.farea[A], mesh.farea[B]) / FACE_TO_STRIP
    exact = 2 * F._tol() + 2 * noise      # (where the mesh's own corners must agree)
    tol = exact + FILLET_TOL
    flat = math.cos(math.radians(1))
    # where both faces run alongside each other, along the edge
    tA, tB = mesh.pts[mesh.fverts[A]] @ d, mesh.pts[mesh.fverts[B]] @ d
    t0, t1 = max(tA.min(), tB.min()), min(tA.max(), tB.max())
    # each face's width across the edge
    width = {X: float(np.ptp(mesh.pts[mesh.fverts[X]] @ np.cross(mesh.fn[X], d))) for X in (A, B)}
    if t1 <= t0:
        return None

    def judge(gs, sign):
        """For facets gs: the radius each one's place suggests if it could be part of
        the band, else nan."""
        gs = np.asarray(gs, dtype=int)
        ok = ~used[gs] & (gs != A) & (gs != B) & (mesh.farea[gs] <= small)
        P, first = corners.of(gs)
        hA, hB = (P - plane_pt[A]) @ n1, (P - plane_pt[B]) @ n2
        absA, absB = np.maximum.reduceat(np.abs(hA), first), np.maximum.reduceat(np.abs(hB), first)
        # (not by its normal: on a mesh with rounded coordinates a sliver lying on the
        # fillet can face any way at all; its corners are what can be trusted. But a
        # piece of either face itself, cut off by the triangulation, isn't band, nor a
        # sliver lying in either face)
        ok &= ~((mesh.fn[gs] @ n1 > flat) & (absA <= tol)) & ~((mesh.fn[gs] @ n2 > flat) & (absB <= tol))
        ok &= (absA > exact) & (absB > exact)
        # it lies between the faces, on the side the ball rolls on
        ok &= (np.maximum.reduceat(sign * hA, first) <= tol) & (np.maximum.reduceat(sign * hB, first) <= tol)
        # and its corners lie on one fillet: each gives the radius of the fillet it would
        # lie on (corners on a wall further off, or a rounded edge along face A's top,
        # give radii far off the fillet's)
        r, bad = _corner_radii(P, p, w_of[sign], d)
        ok &= ~np.maximum.reduceat(bad, first)
        r = np.where(bad, 1.0, r)
        lo, hi = np.minimum.reduceat(r, first), np.maximum.reduceat(r, first)
        ok &= hi <= math.exp(CORNER_AGREE) * np.maximum(lo, 1e-9)
        return np.where(ok, np.add.reduceat(r, first) / corners.count[gs], np.nan)

    def spread(start, accept):
        """The facets reached from start through facets accept() takes (a whole frontier
        at a time)."""
        group = set(int(g) for g in start)
        frontier = list(group)
        while frontier:
            cand = np.array(sorted({h for g in frontier for h in nbrs[g]} - group), dtype=int)
            if not len(cand):
                break
            frontier = cand[accept(cand)].tolist()
            group.update(frontier)
        return np.array(sorted(group), dtype=int)

    w_of = {sg: -sg * (n1 + n2) / (1 + n1 @ n2) for sg in (1.0, -1.0)}
    touched = []
    guess = cache                   # (-1: not judged yet; nan: no)

    def _fillet_side(sign):
        def guesses(gs):
            todo = gs[guess[gs] == -1]
            if len(todo):
                guess[todo] = judge(todo, sign)
                touched.extend(todo.tolist())
            return guess[gs]

        around = np.array(sorted(nbrs[A]), dtype=int)
        seeds = around[np.isfinite(guesses(around))] if len(around) else around
        if not len(seeds):
            return None
        seeds = seeds[np.argsort(guess[seeds], kind="stable")]
        # The facets along face A suggesting about the same radius, grown through those
        # agreeing with them: each such group on its own (beside the fillet, face A
        # can border others, a rounded edge along its top, say), until one reaches B
        # (each group spans no more than DEPTH_AGREE: radii creeping up step by step
        # round a curved edge aren't one fillet)
        groups, first = [], 0
        logs = np.log(guess[seeds])
        for k in range(1, len(seeds) + 1):
            if k == len(seeds) or logs[k] - logs[first] > DEPTH_AGREE:
                groups.append(seeds[first:k])
                first = k
        band = None
        for group in sorted(groups, key=len, reverse=True):
            middle = float(np.median(guess[group]))
            # (the faces must be long and wide enough for a fillet of about this radius:
            # checked now, before growing the band, with the radius only roughly known)
            rough = math.exp(-DEPTH_AGREE) * middle
            if t1 - t0 < FACE_LENGTH * rough or min(width[A], width[B]) < FACE_WIDTH * rough:
                continue
            grown = spread(group, lambda gs: np.abs(np.log(guesses(gs) / middle)) <= DEPTH_AGREE)
            if any(B in nbrs[g] for g in grown):
                band, seeds = grown, group
                break
        if band is None:
            return None
        w = w_of[sign]
        V = mesh.pts[np.unique(np.concatenate([mesh.fverts[g] for g in band]))]
        t = V @ d
        # (a straight fillet is usually cut into strips running its whole length: its
        # corners all lie at its two ends, where the faces end too)
        core = V[(t >= t0 - exact) & (t <= t1 + exact)]
        if len(core) < 4:
            return None
        fit = _radius(core, p, w, d, noise, vouched)
        if fit is None:
            return None
        r, err, allowed = fit
        if t1 - t0 < FACE_LENGTH * r:
            return None                # (faces too short along the edge to tell a fillet by)
        if min(width[A], width[B]) < FACE_WIDTH * r:
            return None                # (a long narrow strip of a finely cut bore is no face)
        a = p + r * w
        # (it may run on past where one face ends, into whatever rounds that end off,
        # by about its radius: strips running its whole length end there; but only
        # facets also reaching into the stretch the faces share)
        run_on = r + exact
        in_band = np.zeros(len(mesh.farea), bool)
        in_band[band] = True

        def on(gs):
            # on the cylinder, and not far past where both faces run (a fillet stops
            # where the edge it rounds off does; whatever carries on round a corner from
            # there is another surface)
            gs = np.asarray(gs, dtype=int)
            P, first = corners.of(gs)
            t = P @ d
            Q = P - a
            Q -= np.outer(Q @ d, d)
            off = np.maximum.reduceat(np.abs(np.linalg.norm(Q, axis=1) - r), first)
            lo, hi = np.minimum.reduceat(t, first), np.maximum.reduceat(t, first)
            # (reaching into the stretch both faces run along: the run-on is the end of a
            # strip running the fillet's length, not a facet of a rounded corner beyond)
            inside = (lo < t1 - exact) & (hi > t0 + exact)
            return in_band[gs] & inside & (lo >= t0 - run_on) & (hi <= t1 + run_on) & (off <= allowed)

        start = seeds[on(seeds)]
        facets = spread(start, on) if len(start) else np.zeros(0, int)
        if not len(facets) or not any(B in nbrs[g] for g in facets):
            return None
        # it must bend on the way (one flat strip with its edges on both faces is a
        # chamfer), and run along its edge for at least about its radius (a narrow slice
        # of a rounded corner or of a rounded edge following a curve fits a cylinder of
        # the same radius too)
        if len(np.unique(np.round(mesh.fn[facets], 2), axis=0)) < 2:
            return None
        t = mesh.pts[np.unique(np.concatenate([mesh.fverts[g] for g in facets]))] @ d
        if np.ptp(t) < FACE_LENGTH * r:
            return None
        # and turn round its axis no further than from one face to the other (more is
        # some other rounded surface it ran on into); nor be part of a bigger cylinder,
        # a hole or pin carrying on round past the faces, through facets on it
        e1 = -w_of[sign] / np.linalg.norm(w_of[sign])
        e2 = np.cross(d, e1)
        turn = math.pi - math.acos(max(-1.0, min(1.0, float(n1 @ n2))))

        def span(gs):
            P = mesh.pts[np.unique(np.concatenate([mesh.fverts[g] for g in gs]))] - a
            ang = np.arctan2(P @ e2, P @ e1)
            # (the smallest arc holding them all, wherever it starts)
            ang = np.sort(np.mod(ang, 2 * math.pi))
            gaps = np.diff(np.r_[ang, ang[0] + 2 * math.pi])
            return 2 * math.pi - gaps.max()

        if span(facets) > turn + math.radians(SPAN_SLACK_DEG):
            return None

        def on_cylinder(gs):
            gs = np.asarray(gs, dtype=int)
            P, first = corners.of(gs)
            Q = P - a
            Q -= np.outer(Q @ d, d)
            off = np.maximum.reduceat(np.abs(np.linalg.norm(Q, axis=1) - r), first)
            return ~used[gs] & (gs != A) & (gs != B) & (off <= allowed)

        whole = spread(facets, on_cylinder)
        if len(whole) > len(facets) and span(whole) > turn + math.radians(SPAN_SLACK_DEG):
            return None
        return facets, sign, a, d, r, err

    for sign in (1.0, -1.0):        # convex edge (ball inside the part) or concave (outside)
        try:
            got = _fillet_side(sign)
        finally:
            cache[touched] = -1.0
            touched.clear()
        if got is not None:
            return got
    return None


def _corner_radii(P, p, w, d):
    """For each point, the radius of the circle touching both faces that passes through
    it on the side facing the faces' corner, as a fillet's points do (the ball's centre
    at p + r w across the edge line; the other circle through it, smaller, has it on the
    far side), and whether there is none (the point lies outside the faces' wedge)."""
    Q = P - p
    Q -= np.outer(Q @ d, d)
    qa, qb, qc = w @ w - 1.0, -2.0 * (Q @ w), np.einsum("ij,ij->i", Q, Q)
    disc = qb * qb - 4 * qa * qc
    # (a corner where the fillet meets a face is on both circles at once, a double root:
    # a hair outside the wedge, from rounding, still counts. Moving a point by e
    # changes the discriminant by about 8 r e |qa|.)
    double = np.abs(qb / (2 * qa))
    root = np.sqrt(np.maximum(disc, 0))
    r = np.maximum((-qb - root) / (2 * qa), (-qb + root) / (2 * qa))
    bad = (disc < -8 * double * abs(qa) * WEDGE_SLACK) | ~(r > 0) | ~np.isfinite(r)
    return r, bad


def _radius(V, p, w, d, noise, vouched):
    """(radius, worst corner gap, gap allowed) of the rolling-ball fillet through corners
    V, or None if they don't lie on one (with vouched radii: on one of those)."""
    Q = V - p
    Q -= np.outer(Q @ d, d)          # position across the edge line
    # A corner q lies on the ball of radius r when |q - r w| = r: a quadratic in r
    qa, qb, qc = w @ w - 1.0, -2.0 * (Q @ w), np.einsum("ij,ij->i", Q, Q)
    if abs(qa) < 1e-12:
        return None
    disc = qb * qb - 4 * qa * qc
    ok = disc >= 0
    if ok.sum() < 3:
        return None
    root = np.sqrt(disc[ok])
    radii = np.r_[(-qb[ok] - root) / (2 * qa), (-qb[ok] + root) / (2 * qa)]
    radii = radii[radii > 1e-3]
    if not len(radii):
        return None

    def gap(r):
        # (most of the corners, not every one: where the band meets a rounded corner or
        # another surface at its ends, a few stray; the fillet is then cut back to the
        # facets lying on it)
        return float(np.percentile(np.abs(np.linalg.norm(Q - r * w, axis=1) - r), ON_SHARE))

    # (each corner gives two radii, the right one shared by most corners: the median of
    # the radii in the densest cluster)
    radii = np.sort(radii)
    cluster = np.split(radii, np.nonzero(np.diff(radii) > max(1e-3, 0.01 * np.median(radii)))[0] + 1)
    r = float(np.median(max(cluster, key=len)))
    # (the band's other facets may stray as far as its middle does, give or take: on a
    # clean mesh that keeps the rounded corners at its ends out, on a rough one it
    # still lets the whole fillet in)
    floor = 2 * F._tol() + 2 * noise

    def cap(limit, radius):
        # (a small fillet is held tighter: a hundredth of a millimetre is a lot of a
        # half-millimetre fillet)
        return min(limit, REL_TOL * radius) + 2 * noise

    err = gap(r)
    # One radius is usually used all over a design: a radius found with confidence
    # elsewhere in the part, close to this band's own, is taken if its corners lie
    # nearly as close to it (or close enough, where they fit no fillet on their own)
    for rv in sorted(vouched, key=lambda x: abs(x - r)):
        if abs(rv - r) > PRIOR_RANGE * rv:
            continue
        e = gap(rv)
        if e <= cap(PRIOR_TOL, rv) and (err > cap(FILLET_TOL, r) or e <= max(2 * err, err + FILLET_TOL / 2)):
            return rv, e, min(cap(max(GROW_TOL, PRIOR_TOL), rv), max(GROW_FACTOR * e, floor))
    if err <= cap(FILLET_TOL, r):
        return r, err, min(cap(GROW_TOL, r), max(GROW_FACTOR * err, floor))
    return None


def _feature(mesh, nbrs, A, B, a, d, r, sign, facets):
    """The fillet's feature (cut to the band's outline), or None if it can't be made."""
    # the faces it joins must be real faces, clearly bigger than its strips (a few
    # facets of a finely cut pin can fit some other cylinder too)
    if min(mesh.farea[A], mesh.farea[B]) < FACE_TO_STRIP * mesh.farea[facets].max():
        return None
    model = Revolved(a, d, line=(r, 0.0))
    F._loose, F._anchored = True, True
    try:
        return F._feature(mesh, model, facets, sign > 0)
    finally:
        F._loose, F._anchored = False, False
