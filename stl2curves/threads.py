"""Screw threads: helical surfaces found and rebuilt exactly.

A thread is swept by its profile (the V of the flanks, with the crest and root between)
turning about an axis while it moves along it. Every surface normal of such a sweep
lies in one "linear line complex" (Pottmann and Randrup, "Rotational and helical surface
approximation for reverse engineering", 1998): a normal n at a point x satisfies

    c . (x cross n) + cbar . n = 0

for the six numbers (c, cbar) of the motion, so the axis (direction c, through
c cross cbar) and the pitch (c . cbar per radian) come from one small eigenvalue
problem over the facets' normals, finished by least squares on the angles.

The normals of a finely cut, rounded-off mesh are only good to a degree or so, which
fixes the axis roughly but hardly the pitch. The corners are good to the file's
precision, and a thread's corners are special: its crest and root edges run round
every turn, so the corners bunch up at a few radii (fixing the axis exactly: each
bunch is a cylinder round it), and within a bunch they lie on a few helix lines (the
pitch is the one that lines them up). Turned back along the helix into one plane
through the axis, the corners of every turn then land on the same few points: the
corners of the profile.

Any smooth surface agrees with some screw motion over a small enough patch, so a
thread must also wind most of the way round its axis (a finely meshed bend agrees with
one over a narrow arc) and be no deeper than its lead.

The thread is rebuilt one face per straight piece of its profile (each flank, the
crest, the root): the piece swept along the helix, a ruled B-spline surface between
two helices (exact for a straight piece), cut to the outline of its facets.
"""

import math

import numpy as np

MIN_FACETS = 30        # a thread has at least this many facets
FIT_DEG = 2.0          # a seed's facet normals this close to a screw motion ...
SEED_SHARE = 0.25      # ... this many of them: the motion is worth a closer look
AGREE_DEG = 3.0        # facet normals this close to the motion may belong to the thread
SMALLEST_SEED = 1.0    # mm: radius of the smallest balls of facets tried as seeds
MAX_FITS = 300         # seeds tried at most, per part
BAND = 0.02            # mm: radii binned this finely to find the crest's and root's corners
MIN_GATHER = 0.8       # share of a crest's or root's corners on a few helix lines
ON_TOL = 0.01          # mm (plus the file's rounding): corners this near lie on the thread
MIN_TURNS = 1.0        # the facets must wind round at least this far
PROFILE_TOL = 0.004    # mm: profile simplified to straight pieces within this
MIN_AROUND = 270       # deg: the facets must reach at least this far round the axis
MAX_FACET_TURN = 90    # deg: a facet of a piece spans less than this along its helix (a plain
                       # bore a whole number of leads long has every corner on the thread)


class Helical:
    """Profile (z, rho) pairs over one lead, swept by the screw motion: axis through a
    along unit d, rising `pitch` per radian. Points outside the material are positive."""

    def __init__(self, a, d, pitch, profile, inside):
        self.d = d / np.linalg.norm(d)
        self.e1 = np.cross(self.d, [1.0, 0, 0] if abs(self.d[0]) < 0.9 else [0, 1.0, 0])
        self.e1 /= np.linalg.norm(self.e1)
        self.e2 = np.cross(self.d, self.e1)
        self.a = np.asarray(a, float) - (np.asarray(a, float) @ self.d) * self.d
        self.pitch = float(pitch)
        self.lead = 2 * math.pi * abs(self.pitch)
        self.profile = np.asarray(profile, float)       # closed over one lead: z from 0 to lead
        self.inside = inside                            # material nearer the axis (a bolt)
        self.kind = "thread"

    def local(self, p):
        """(rho, z reduced to the profile's plane, angle) of points p."""
        v = np.atleast_2d(p) - self.a
        z = v @ self.d
        w = v - np.outer(z, self.d)
        u = np.arctan2(w @ self.e2, w @ self.e1)
        return np.linalg.norm(w, axis=1), np.mod(z - self.pitch * u, self.lead), u

    def signed(self, p):
        rho, z0, _ = self.local(p)
        return self._profile_distance(z0, rho)

    def _profile_distance(self, z0, rho):
        P = self.profile
        best = np.full(len(z0), np.inf)
        sign = np.ones(len(z0))
        for shift in (-self.lead, 0.0, self.lead):
            A, B = P[:-1] + [shift, 0], P[1:] + [shift, 0]
            q = np.c_[z0, rho][:, None]
            ab = B - A
            t = np.clip(((q - A) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-300), 0, 1)
            foot = A + ab * t[..., None]
            d = np.linalg.norm(q - foot, axis=-1)
            k = d.argmin(axis=1)
            dk = d[np.arange(len(z0)), k]
            # side: rho above the profile at that z is away from a bolt's material
            above = rho - foot[np.arange(len(z0)), k, 1]
            better = dk < best
            best[better] = dk[better]
            sign[better] = np.where(above[better] >= 0, 1.0, -1.0)
        return best * (sign if self.inside else -sign)

    def point(self, z0, rho, u):
        """The point of the surface over profile point (z0, rho) at angle u."""
        z = z0 + self.pitch * u
        return self.a + z * self.d + rho * (math.cos(u) * self.e1 + math.sin(u) * self.e2)


class Flank:
    """One straight piece of a thread's profile, from (z, rho) p to q, swept along the
    thread's helix: a flank (a helicoid), or the crest or root between (a helical band
    round a cylinder). A point of it is (s, U): s mm along the piece from p, on the
    helix turned U radians round (unwrapped: U tells the turns apart). Points away
    from the material are positive."""

    kind = "thread"

    def __init__(self, thread, p, q):
        self.thread = thread
        self.p, self.q = np.asarray(p, float), np.asarray(q, float)
        self.length = float(np.linalg.norm(self.q - self.p))
        self.t = (self.q - self.p) / self.length            # along the piece, as (z, rho)
        side = 1.0 if thread.inside else -1.0
        # (the profile runs with z rising, the material below it for a bolt)
        self.n = side * np.array([-self.t[1], self.t[0]])   # away from the material
        self.a, self.d = thread.a, thread.d

    def place(self, x):
        """(s, distance off the piece's line in the profile's plane, U) of points x,
        turned back along the helix onto the piece (the nearest of a lead either way)."""
        h = self.thread
        v = np.atleast_2d(x) - h.a
        z = v @ h.d
        w = v - np.outer(z, h.d)
        rho = np.linalg.norm(w, axis=1)
        z0 = np.mod(z - h.pitch * np.arctan2(w @ h.e2, w @ h.e1), h.lead)
        best = None
        for shift in (-h.lead, 0.0, h.lead):
            q = np.c_[z0 + shift, rho] - self.p
            s, off = q @ self.t, q @ self.n
            gap = np.hypot(np.maximum(0, np.maximum(-s, s - self.length)), off)  # to the piece
            here = (s, off, z0 + shift, gap)
            if best is None:
                best = here
            else:
                better = gap < best[3]
                best = tuple(np.where(better, x1, x0) for x0, x1 in zip(best, here))
        s, off, zp, _ = best
        return s, off, (z - zp) / h.pitch

    def point(self, s, U):
        """Points of the surface at (s, U)."""
        h = self.thread
        s, U = np.broadcast_arrays(np.asarray(s, float), np.asarray(U, float))
        z = self.p[0] + s * self.t[0] + h.pitch * U
        rho = self.p[1] + s * self.t[1]
        r = np.cos(U)[..., None] * h.e1 + np.sin(U)[..., None] * h.e2
        return h.a + z[..., None] * h.d + rho[..., None] * r

    def _normals(self, s, U):
        """Unit normals at (s, U) (away from the material), and the profile plane's own
        normal there (the distance measured in that plane is a hair more than the
        true one: the helix leans it)."""
        h = self.thread
        rho = self.p[1] + s * self.t[1]
        r = np.cos(U)[:, None] * h.e1 + np.sin(U)[:, None] * h.e2
        across = -np.sin(U)[:, None] * h.e1 + np.cos(U)[:, None] * h.e2
        n = np.cross(self.t[0] * h.d + self.t[1] * r, h.pitch * h.d + rho[:, None] * across)
        n /= np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None]
        flat = self.n[0] * h.d + self.n[1] * r
        n[np.einsum("ij,ij->i", n, flat) < 0] *= -1
        return n, flat

    def signed(self, x):
        s, off, U = self.place(x)
        n, flat = self._normals(s, U)
        return off * np.einsum("ij,ij->i", n, flat)

    def normal(self, x):
        s, _, U = self.place(x)
        return self._normals(s, U)[0]


def screw_fit(X, N, w):
    """Screw motion fitted to points X with normals N (weights w): (axis direction c,
    point on the axis, pitch per radian)."""
    o = (X * w[:, None]).sum(0) / w.sum()
    M = np.c_[np.cross(X - o, N), N]
    A = (M * w[:, None]).T @ M
    A11, A12, A22 = A[:3, :3], A[:3, 3:], A[3:, 3:]
    A22i = np.linalg.pinv(A22)
    lam, V = np.linalg.eigh(A11 - A12 @ A22i @ A12.T)
    c = V[:, 0]
    cbar = -A22i @ A12.T @ c
    return c, np.cross(c, cbar) + o, float(c @ cbar)


def misfit_deg(X, N, c, a, pitch):
    """Angle between each normal and the plane square to the screw motion's path there."""
    v = np.cross(c, X - a) + pitch * c
    s = np.abs(np.einsum("ij,ij->i", v, N)) / np.maximum(np.linalg.norm(v, axis=1), 1e-12)
    return np.degrees(np.arcsin(np.clip(s, 0, 1)))


def find(mesh, candidates):
    """[(Helical, facets)] for the threads among the candidate facets. Each connected
    patch of them is tried whole, then (a thread often joins smoothly onto other
    unexplained facets: fillets round the boss, run-outs) in balls of facets spread
    over it, halving in size, until too few facets are left."""
    found, tested, fits = [], [], 0
    for patch in _patches(mesh, np.asarray(candidates)):
        pool = patch
        for seed in _seeds(mesh, patch):
            seed = np.intersect1d(seed, pool)
            if len(seed) < MIN_FACETS:
                continue
            fits += 1
            if fits > MAX_FITS:
                return found
            motion = _motion(mesh, seed)
            if motion is None or any(_same(motion, seed, t) for t in tested):
                continue
            thread, agree = _thread(mesh, pool, seed, *motion)
            if thread is None:
                tested.append((motion, agree))
                continue
            found.append(thread)
            pool = np.setdiff1d(pool, thread[1])
            if len(pool) < MIN_FACETS:
                break
    return found


METRIC_PITCHES = (0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 1.0, 1.25, 1.5,
                  1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)


def name(thread):
    """The thread in plain words, e.g. "M12 x 1.5, right-hand, external"."""
    major = 2 * float(thread.profile[:, 1].max())
    lead = thread.lead
    tpi = 25.4 / lead
    nearest = min(METRIC_PITCHES, key=lambda p: abs(p - lead))
    if abs(nearest - lead) <= 0.01 * lead and abs(major - round(major)) <= 0.4:
        size = f"M{round(major)} x {nearest:g}"
    elif abs(tpi - round(tpi)) <= 0.01 * tpi:
        size = f'{major / 25.4:.3f}" x {round(tpi)} TPI'
    else:
        size = f"dia {major:.2f} x pitch {lead:.3f}"
    hand = "right-hand" if thread.pitch > 0 else "left-hand"
    return f"{size}, {hand}, {'external' if thread.inside else 'internal'}"


def patches(mesh, thread, facets):
    """The thread's facets as features, one per connected stretch of each straight piece
    of its profile (a facet must lie wholly on one piece; any that don't are left as
    they are)."""
    tol = ON_TOL + mesh.noise
    P = thread.profile
    pieces = [Flank(thread, P[i], P[i + 1]) for i in range(len(P) - 1)]
    corners = np.unique(np.concatenate([mesh.fverts[f] for f in facets]))
    on = np.zeros((len(pieces), len(mesh.pts)), bool)
    U = np.zeros((len(pieces), len(mesh.pts)))
    for k, piece in enumerate(pieces):
        s, off, U[k, corners] = piece.place(mesh.pts[corners])
        on[k, corners] = (np.abs(off) <= tol) & (s >= -tol) & (s <= piece.length + tol)
    groups = [[] for _ in pieces]
    for f in facets:
        homes = [k for k in range(len(pieces)) if on[k, mesh.fverts[f]].all()
                 and np.ptp(U[k, mesh.fverts[f]]) < math.radians(MAX_FACET_TURN)]
        if len(homes) == 1:
            groups[homes[0]].append(int(f))
    out = []
    for piece, group in zip(pieces, groups):
        # each connected stretch of a piece is a face of its own (other features may
        # have taken facets out of the middle, and the thread's run-outs leave scraps)
        left = set(group)
        while left:
            part = _reach(mesh, [min(left)], list(left))
            left -= set(part.tolist())
            if len(part) >= 3:
                out.append(_piece_feature(mesh, thread, piece, part))
    return out


def _piece_feature(mesh, thread, piece, fids):
    from .features import Feature
    tids = np.concatenate([mesh.ftris[f] for f in fids])
    T = mesh.tris[tids]
    # how far each facet sits from the true surface, at its edge midpoints
    corners = mesh.pts[T]
    gap = piece.signed(((corners + corners[:, [1, 2, 0]]) / 2).reshape(-1, 3)).reshape(-1, 3)
    area = mesh.tarea[tids]
    s, _, U = piece.place(mesh.pts[np.unique(T)])
    rho = thread.profile[:, 1]
    if abs(piece.t[1]) > 0.1:
        part = "flank"
    elif piece.p[1] > (rho.min() + rho.max()) / 2:
        part = "crest" if thread.inside else "root"
    else:
        part = "root" if thread.inside else "crest"
    return Feature(piece, "screw thread", f"{name(thread)}: {part}", convex=thread.inside,
                   kind="trimmed", facets=fids,
                   change=-float((area * gap.mean(axis=1)).sum()),
                   tolerance=0.5 * float((area * np.abs(gap).mean(axis=1)).sum()) + 1e-3,
                   worst=float(np.abs(gap).max()),
                   u0=float(U.min()), span=float(np.ptp(U)), lo=float(s.min()), hi=float(s.max()))


def _seeds(mesh, patch):
    """The patch itself, then balls of its facets round points spread over it, the
    balls halving in radius down to SMALLEST_SEED."""
    yield patch
    X = mesh.fcent[patch]
    r = 0.25 * float(np.ptp(X, axis=0).max())
    while r >= SMALLEST_SEED:
        for i in _spread(X, r):
            yield patch[np.linalg.norm(X - X[i], axis=1) <= r]
        r /= 2


def _spread(X, r):
    """Points of X about r apart covering them all, farthest first (lazily)."""
    i = 0
    dist = np.full(len(X), np.inf)
    while True:
        yield i
        dist = np.minimum(dist, np.linalg.norm(X - X[i], axis=1))
        i = int(dist.argmax())
        if dist[i] <= r:
            return


def _same(motion, seed, tested):
    """Whether a motion (and seed) was already tested: the same axis, and the seed
    mostly among the facets tried with it."""
    (c, a, _), ((c2, a2, _), facets) = motion, tested
    if abs(c @ c2) < math.cos(math.radians(1)):
        return False
    off = a - a2
    if np.linalg.norm(off - (off @ c2) * c2) > 0.1:
        return False
    return len(np.intersect1d(seed, facets)) >= 0.5 * len(seed)


def _motion(mesh, facets):
    """(axis direction, point on the axis, pitch per radian) of the screw motion most
    of the facets' normals agree with, or None. The eigenvalue fit gives the start,
    outliers dropped; least squares on the angles themselves finishes it (the
    eigenvalue fit's own measure favours an axis near the facets, where the motion
    barely moves them, so on part of a thread its axis wanders off)."""
    X, N, w = mesh.fcent[facets], mesh.fn[facets], mesh.farea[facets]
    keep = np.ones(len(facets), bool)
    for _ in range(10):
        if keep.sum() < MIN_FACETS:
            return None
        c, a, pitch = screw_fit(X[keep], N[keep], w[keep])
        err = misfit_deg(X, N, c, a, pitch)
        new = err < max(FIT_DEG, 3 * float(np.median(err[keep])))
        if (new == keep).all():
            break
        keep = new
    c, a, pitch = _angle_fit(X, N, w, c, a, pitch)
    if np.mean(misfit_deg(X, N, c, a, pitch) < FIT_DEG) < SEED_SHARE:
        return None
    if np.median(_radius(X, c, a)) > 2 * np.ptp(X, axis=0).max():
        return None                     # an axis far off: the facets don't wind round it
    return c, a, pitch


def _angle_fit(X, N, w, c, a, pitch):
    """Screw motion refined to the angles between the normals and the motion's paths
    (robust least squares: a degree or two of disagreement is the facets' own)."""
    from scipy.optimize import least_squares
    d0, e1, e2 = _frame(c)
    scale = np.sqrt(w / w.mean())

    def unpack(x):
        d = d0 + x[0] * e1 + x[1] * e2
        return d / np.linalg.norm(d), a + x[2] * e1 + x[3] * e2, pitch + x[4]

    def residual(x):
        d, p, k = unpack(x)
        v = np.cross(d, X - p) + k * d
        s = np.einsum("ij,ij->i", v, N) / np.maximum(np.linalg.norm(v, axis=1), 1e-12)
        return np.degrees(np.arcsin(np.clip(s, -1, 1))) * scale

    fit = least_squares(residual, np.zeros(5), loss="soft_l1", f_scale=1.0, x_scale="jac", max_nfev=100)
    return unpack(fit.x)


def _coaxial(P, c, a):
    """Axis refitted to the corners on the crest and root edges, each set on a cylinder
    round it: the corners are good to the file's precision, where the normals leave
    the axis a degree or two out."""
    from scipy.optimize import least_squares
    for width in (0.15, 0.06, 0.03):    # the bands sharpen as the axis comes right
        rho = _radius(P, c, a)
        radii = _bands(rho, 2, width / 1.5)
        if not radii:
            return c, a
        which = np.full(len(P), -1)
        for k, r in enumerate(radii):
            which[np.abs(rho - r) <= width] = k
        sel = which >= 0
        if sel.sum() < 10:
            return c, a
        d0, e1, e2 = _frame(c)

        def unpack(x):
            d = d0 + x[0] * e1 + x[1] * e2
            return d / np.linalg.norm(d), a + x[2] * e1 + x[3] * e2

        def residual(x):
            d, p = unpack(x)
            return _radius(P[sel], d, p) - x[4:][which[sel]]

        fit = least_squares(residual, np.r_[0, 0, 0, 0, radii], loss="soft_l1", f_scale=0.01, x_scale="jac")
        c, a = unpack(fit.x)
    return c, a


def _bands(rho, most, cell=BAND):
    """Radii the points bunch up at (a thread's crest and root edges), densest first."""
    edges = np.arange(rho.min(), rho.max() + cell, cell)
    if len(edges) < 2:
        return []
    count = np.histogram(rho, edges)[0]
    out = []
    for i in np.argsort(count)[::-1]:
        if count[i] < max(10, 0.2 * count.max()) or len(out) == most:
            break
        r = edges[i] + cell / 2
        if all(abs(r - s) > 3 * cell for s in out):
            out.append(r)
    return out


def _radius(X, c, a):
    """Distances of points X from the axis through a along unit c."""
    v = X - a
    return np.linalg.norm(v - np.outer(v @ c, c), axis=1)


def _frame(d):
    d = np.asarray(d, float) / np.linalg.norm(d)
    e1 = np.cross(d, [1.0, 0, 0] if abs(d[0]) < 0.9 else [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    return d, e1, np.cross(d, e1)


def _thread(mesh, pool, seed, c, a, pitch):
    """(thread, facets) about this axis through the seed, and the facets tried. The
    normals only give the axis (and the pitch roughly): the pitch and the profile come
    from the corners, which lie on the thread to the file's precision."""
    X, N, W = mesh.fcent[pool], mesh.fn[pool], mesh.farea[pool]
    reach = 1.5 * _radius(mesh.fcent[seed], c, a).max()
    agree = seed
    for _ in range(6):                  # facets agreeing, the motion refitted to them, ...
        near = (misfit_deg(X, N, c, a, pitch) < AGREE_DEG) & (_radius(X, c, a) <= reach)
        if near.sum() < MIN_FACETS:
            return None, agree
        if _around(X[near], c, a) < MIN_AROUND:
            # (a thread winds round its axis: a smooth bend's facets can agree with
            # some screw motion, but only over a narrow arc of it)
            return None, pool[near]
        if np.array_equal(pool[near], agree):
            break
        agree = pool[near]
        c, a, pitch = _angle_fit(X[near], N[near], W[near], c, a, pitch)
        # (fitted to only the facets agreeing so far, the axis can settle a degree or
        # so out and the pitch far out, keeping the rest out; the corners pull them back)
        c, a = _coaxial(_corners(mesh, agree), c, a)
        pitch = _corner_pitch(_corners(mesh, agree), c, a) or pitch
    P = _corners(mesh, agree)
    pitch = _corner_pitch(P, c, a)
    if pitch is None:
        return None, agree
    model = Helical(a, c, pitch, [[0, 1], [1, 1]], True)
    prof = _profile(model, P)
    if prof is None:
        return None, agree
    for near in (0.05, 0.02):           # the axis from the normals is a little off at first
        rho, z0, _ = model.local(P)
        trial = Helical(model.a, model.d, model.pitch, prof, True)
        on = np.abs(trial._profile_distance(z0, rho)) <= near
        if on.sum() < 3 * len(prof):
            return None, agree
        model, prof = _refine(model, prof, P[on])
    # grown onto every facet lying on it, refitted as they join (fitted only to the
    # facets whose normals agreed, it can be a hair out over a long thread)
    for near in (0.05, 0.02):
        trial = Helical(model.a, model.d, model.pitch, prof, True)
        grown = _lying_on(mesh, pool, agree, trial, near)
        if len(grown) < MIN_FACETS:
            return None, agree
        model, prof = _refine(model, prof, _corners(mesh, grown))
    if (np.diff(prof[:, 0]) <= 0).any():
        return None, agree              # the profile folded over itself: not a thread
    if np.ptp(prof[:, 1]) > model.lead:
        return None, agree              # deeper than its lead: no thread is
    prof = _simplify(prof, PROFILE_TOL)
    # a crest or root flat to within the file's rounding: exactly flat (only then: a
    # thread's own crest can slope or bow by a few microns, and flattening it would
    # push its corners off)
    ring = prof[:-1].copy()
    for i in range(len(ring)):
        j = (i + 1) % len(ring)
        if abs(ring[j, 1] - ring[i, 1]) <= max(2 * mesh.noise, 0.001):
            ring[i, 1] = ring[j, 1] = (ring[i, 1] + ring[j, 1]) / 2
    prof = np.r_[ring, ring[:1] + [model.lead, 0]]
    model = Helical(model.a, model.d, model.pitch, prof, True)
    facets = _lying_on(mesh, pool, agree, model, ON_TOL + mesh.noise)
    if len(facets) < MIN_FACETS:
        return None, agree
    Q = mesh.pts[np.unique(np.concatenate([mesh.fverts[f] for f in facets]))]
    if np.ptp((Q - model.a) @ model.d) < MIN_TURNS * model.lead or _around(Q, model.d, model.a) < MIN_AROUND:
        return None, agree              # the facets must wind round the axis
    # which side is material: facets' normals point away from it
    X, N = mesh.fcent[facets], mesh.fn[facets]
    radial = X - model.a - np.outer((X - model.a) @ model.d, model.d)
    inside = (np.einsum("ij,ij->i", radial, N) * mesh.farea[facets]).sum() > 0
    return (Helical(model.a, model.d, model.pitch, prof, inside), facets), agree


def _around(X, c, a):
    """How far (degrees) the points reach round the axis through a along c: all the way
    less the widest angle between them."""
    _, e1, e2 = _frame(c)
    v = X - a
    u = np.sort(np.arctan2(v @ e2, v @ e1))
    return 360.0 - math.degrees(np.diff(np.r_[u, u[0] + 2 * math.pi]).max())


def _corners(mesh, facets):
    return mesh.pts[np.unique(np.concatenate([mesh.fverts[f] for f in facets]))]


def _lying_on(mesh, pool, start, model, tol):
    """Facets of the pool with every corner within tol of the model, reached from those
    of `start` through shared edges."""
    corners = np.unique(np.concatenate([mesh.fverts[f] for f in pool]))
    off = np.full(len(mesh.pts), np.inf)
    off[corners] = np.abs(model.signed(mesh.pts[corners]))
    return _reach(mesh, start, [f for f in pool if (off[mesh.fverts[f]] <= tol).all()])


def _reach(mesh, start, allowed):
    """Facets of `allowed` reached from those of `start` among them through shared edges."""
    allowed = set(np.asarray(allowed).tolist())
    todo = [f for f in np.asarray(start).tolist() if f in allowed]
    seen = set(todo)
    while todo:
        f = todo.pop()
        for g in _adjacent(mesh)[f]:
            if g in allowed and g not in seen:
                seen.add(g)
                todo.append(g)
    return np.array(sorted(seen), int)


def _corner_pitch(P, c, a):
    """_pitch of corners P about the axis through a along c."""
    d, e1, e2 = _frame(c)
    v = P - a
    z = v @ d
    w = v - np.outer(z, d)
    return _pitch(z, np.arctan2(w @ e2, w @ e1), np.linalg.norm(w, axis=1))


def _pitch(z, u, rho):
    """Pitch per radian (signed: right- or left-handed) of the helix lines the corners
    run along, or None. The crest and root edges of a thread run round every turn, so
    the corners bunch up at a few radii; in each bunch, turned back along the right
    helix, they fall on a few lines (one per edge)."""
    best = (0.0, None)
    for r in _bands(rho, 3):
        band = np.abs(rho - r) <= 2 * BAND
        best = max(best, _lines(z[band], u[band], rho.max()), key=lambda b: b[0])
    return best[1] if best[0] >= MIN_GATHER else None


def _lines(z, u, radius):
    """(share of the points on at most a few helix lines, pitch per radian) for the
    best lead from 0.1 mm to half the radius, either hand: turned back along the
    right helix, each line's points share one phase."""
    turn = u / (2 * math.pi)

    def scan(leads, hands, bins, top):
        best = (0.0, None)
        for hand in hands:
            for chunk in np.array_split(leads, max(1, len(leads) // 200)):
                phase = np.mod(z[None, :] / chunk[:, None] - hand * turn[None, :], 1.0)
                cells = np.minimum((phase * bins).astype(int), bins - 1)
                cells += bins * np.arange(len(chunk))[:, None]
                hist = np.bincount(cells.ravel(), minlength=bins * len(chunk)).reshape(len(chunk), bins)
                share = np.sort(hist, axis=1)[:, -top:].sum(axis=1) / len(z)
                k = int(np.argmax(share))
                if share[k] > best[0]:
                    best = (float(share[k]), hand * float(chunk[k]))
        return best

    leads = np.linspace(0.1, max(0.5 * radius, 0.2), 4000)
    share, lead = scan(leads, (1, -1), 32, 4)
    if lead is None:
        return 0.0, None
    step = leads[1] - leads[0]
    fine = np.linspace(abs(lead) - step, abs(lead) + step, 101)
    _, finer = scan(fine, (1 if lead > 0 else -1,), 128, 8)
    return share, (finer if finer is not None else lead) / (2 * math.pi)


def _profile(model, P):
    """The profile's corners (z along one lead, rho), closed, or None. Turned back
    along the helix into one plane through the axis, the corners gather where the
    mesh's helix lines run (the crest and root edges, every turn putting its corners
    in the same places). The gatherings are taken down to the biggest drop in size:
    past it are the scattered corners of run-outs and of other surfaces."""
    rho, z0, _ = model.local(P)
    cell = max(0.005, model.lead / 150)
    cols = int(math.ceil(model.lead / cell))
    zi = np.minimum((z0 / cell).astype(int), cols - 1)
    ri = np.floor(rho / cell).astype(int)
    cells, inv = np.unique(zi * 100000 + ri, return_inverse=True)
    inv = inv.ravel()
    at = {k: i for i, k in enumerate(cells.tolist())}
    parent = list(range(len(cells)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, k in enumerate(cells.tolist()):  # join touching cells (z wraps round)
        cz, cr = divmod(k, 100000)
        for dz in (-1, 0, 1):
            for dr in (-1, 0, 1):
                j = at.get(((cz + dz) % cols) * 100000 + cr + dr)
                if j is not None:
                    parent[root(i)] = root(j)
    group = np.array([root(i) for i in range(len(cells))])[inv]
    ids, size = np.unique(group, return_counts=True)
    order = np.argsort(size)[::-1][:40]
    size = size[order]
    if len(size) < 3:
        return None
    drop = size[1:-1] / size[2:]        # cutting after the 2nd, 3rd, ... gathering
    keep = int(np.argmax(drop)) + 2
    if drop[keep - 2] < 3:
        return None                     # no clear cut between edges and scatter
    prof = []
    for g in ids[order[:keep]]:
        sel = group == g
        angle = 2 * math.pi * z0[sel] / model.lead            # (a circular mean: may straddle 0)
        zc = math.atan2(np.sin(angle).mean(), np.cos(angle).mean()) % (2 * math.pi) * model.lead / (2 * math.pi)
        prof.append((zc, float(np.median(rho[sel]))))
    prof = np.array(sorted(prof))
    return np.r_[prof, prof[:1] + [model.lead, 0]]          # closed over one lead


def _adjacent(mesh):
    """Facets sharing an edge with each facet, however sharply they meet (cached)."""
    if "all_nbrs" not in mesh.__dict__:
        nbrs = [set() for _ in range(len(mesh.fn))]
        for sides in mesh.edge_facets.values():
            for f in sides:
                nbrs[f].update(int(g) for g in sides if g != f)
        mesh.all_nbrs = [sorted(x) for x in nbrs]
    return mesh.all_nbrs


def _patches(mesh, facets):
    """Groups of the facets joined through shared edges (however sharply they bend)."""
    inside = set(facets.tolist())
    parent = {f: f for f in inside}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for sides in mesh.edge_facets.values():
        mine = [int(f) for f in sides if int(f) in inside]
        for f in mine[1:]:
            a, b = root(mine[0]), root(f)
            if a != b:
                parent[a] = b
    groups = {}
    for f in inside:
        groups.setdefault(root(f), []).append(f)
    return [np.array(sorted(g)) for g in groups.values() if len(g) >= MIN_FACETS]


def _refine(model, prof, P):
    """Axis (tilt and position), pitch and profile corners fitted together (least
    squares) to corners P that lie on the thread."""
    from scipy.optimize import least_squares
    d0, e1, e2 = model.d, model.e1, model.e2
    k = len(prof) - 1                   # corners (the last repeats the first a lead on)

    def build(x):
        d = d0 + x[0] * e1 + x[1] * e2
        a = model.a + x[2] * e1 + x[3] * e2
        pitch = model.pitch * (1 + x[4])
        corners = prof[:-1] + x[5:].reshape(k, 2)
        lead = 2 * math.pi * abs(pitch)
        closed = np.r_[corners, corners[:1] + [lead, 0]]
        return Helical(a, d, pitch, closed, True)

    def residual(x):
        m = build(x)
        rho, z0, _ = m.local(P)
        return m._profile_distance(z0, rho)

    fit = least_squares(residual, np.zeros(5 + 2 * k), x_scale="jac", max_nfev=200)
    m = build(fit.x)
    return m, m.profile


def _simplify(P, tol):
    """Douglas-Peucker on an open polyline (keeps both ends)."""
    keep = np.zeros(len(P), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(P) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = P[i], P[j]
        ab = b - a
        t = np.clip(((P[i + 1:j] - a) @ ab) / max(ab @ ab, 1e-300), 0, 1)
        d = np.linalg.norm(P[i + 1:j] - (a + np.outer(t, ab)), axis=1)
        k = int(np.argmax(d))
        if d[k] > tol:
            m = i + 1 + k
            keep[m] = True
            stack += [(i, m), (m, j)]
    return P[keep]
