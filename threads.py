"""Screw threads: helical surfaces found and rebuilt exactly.

A thread is swept by its profile (the V of the flanks, with the crest and root between)
turning about an axis while it moves along it. Every surface normal of such a sweep
lies in one "linear line complex" (Pottmann and Randrup, "Rotational and helical surface
approximation for reverse engineering", 1998): a normal n at a point x satisfies

    c . (x cross n) + cbar . n = 0

for the six numbers (c, cbar) of the motion, so the axis (direction c, through
c cross cbar) and the pitch (c . cbar per radian) come from one small eigenvalue
problem over the facets' normals. The facets no simple surface explained are fitted
this way, outliers dropped and the fit repeated; a thread is accepted when many facets
agree with one screw motion of real pitch to a fraction of a degree. Its profile is
then read off by turning every corner back along the helix into one plane through the
axis, where the corners of all the turns fall on one curve.

The thread is rebuilt as one surface: the profile swept along the helix (a B-spline
surface, straight between the profile's corners), cut to the outline of its facets.
"""

import math

import numpy as np

MIN_FACETS = 60        # a thread has at least this many facets
MAX_ANGLE_DEG = 0.3    # median disagreement between facet normals and the screw motion
MIN_TURNS = 1.0        # the facets must wind round at least this far
PROFILE_TOL = 0.004    # mm: profile simplified to straight pieces within this


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
    """[(Helical, facets)] for the threads among the candidate facets: each connected
    patch of them is tried on its own (a part may have several threads, and plenty of
    other unexplained facets)."""
    found = []
    for patch in _patches(mesh, np.asarray(candidates)):
        cand = patch
        for _ in range(3):
            if len(cand) < MIN_FACETS:
                break
            thread = _one(mesh, cand)
            if thread is None:
                break
            found.append(thread)
            cand = np.setdiff1d(cand, thread[1])
    return found


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


def _one(mesh, cand):
    X, N, w = mesh.fcent[cand], mesh.fn[cand], mesh.farea[cand]
    keep = np.ones(len(cand), bool)
    for _ in range(10):
        if keep.sum() < MIN_FACETS:
            return None
        c, a, pitch = screw_fit(X[keep], N[keep], w[keep])
        err = misfit_deg(X, N, c, a, pitch)
        new = err < max(0.5, 3 * float(np.median(err[keep])))
        if (new == keep).all():
            break
        keep = new
    if keep.sum() < MIN_FACETS or np.median(err[keep]) > MAX_ANGLE_DEG:
        return None
    lead = 2 * math.pi * abs(pitch)
    rho = np.linalg.norm(np.cross(X[keep] - a, c), axis=1)
    if lead < 0.05 or lead > 0.5 * rho.max():
        return None                     # no real pitch (a surface of revolution) or too steep
    facets = cand[keep]
    # The normals give the axis well but the pitch only roughly (coarse facets): the
    # pitch is the one that lines the crest corners of every turn up in one plane
    model = Helical(a, c, pitch, [[0, 1], [1, 1]], True)
    corners = np.unique(np.concatenate([mesh.fverts[f] for f in facets]))
    v = mesh.pts[corners] - model.a
    z = v @ model.d
    w = v - np.outer(z, model.d)
    rho = np.linalg.norm(w, axis=1)
    u = np.arctan2(w @ model.e2, w @ model.e1)
    crest = rho >= rho.min() + 0.66 * np.ptp(rho)
    if crest.sum() < 8:
        return None
    pitch, coherence = _best_pitch(z[crest], u[crest], rho.max())
    if pitch is None or coherence < 0.6:
        return None
    model = Helical(a, c, pitch, [[0, 1], [1, 1]], True)
    if np.ptp(z) / model.lead < MIN_TURNS:
        return None                     # the facets must wind round the axis at least once
    return _with_profile(mesh, model, facets)


def _best_pitch(z, u, radius):
    """Pitch per radian (signed: right- or left-handed) that brings the points' phases
    z / pitch - angle into line, and how well (1: perfectly)."""
    best = (0.0, None)
    for hand in (1, -1):
        for lead in np.linspace(0.1, 0.5 * radius, 4000):
            phase = 2 * math.pi * z / lead - hand * u
            r = abs(np.exp(1j * phase).mean())
            if r > best[0]:
                best = (r, hand * lead)
    if best[1] is None:
        return None, 0.0
    # refine around the best
    r0, lead0 = best
    hand = 1 if lead0 > 0 else -1
    step = 0.5 * radius / 4000
    for lead in np.linspace(abs(lead0) - step, abs(lead0) + step, 201):
        phase = 2 * math.pi * z / lead - hand * u
        r = abs(np.exp(1j * phase).mean())
        if r > best[0]:
            best = (r, hand * lead)
    return best[1] / (2 * math.pi), best[0]


def _with_profile(mesh, model, facets):
    """The thread's profile. Turned back along the helix into one plane through the
    axis, the corners of every turn gather on the profile's own corners (the crest and
    root edges run along the helix, so each turn puts its corners in the same places);
    the dense gatherings, in order along the axis, are the profile. (Corners where the
    thread runs out fall elsewhere, but only a few in any one place.)"""
    corners = np.unique(np.concatenate([mesh.fverts[f] for f in facets]))
    rho, z0, _ = model.local(mesh.pts[corners])
    cell = max(0.005, model.lead / 150)
    cols = int(math.ceil(model.lead / cell))
    zi = np.minimum((z0 / cell).astype(int), cols - 1)
    ri = np.floor(rho / cell).astype(int)
    key = zi * 100000 + ri
    cells, inv, count = np.unique(key, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    dense = count >= max(3, 0.08 * count.max())
    if dense.sum() < 3:
        return None
    # gatherings: dense cells next to each other (the axial direction wraps round)
    dz, dr = cells[dense] // 100000, cells[dense] % 100000
    parent = list(range(len(dz)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(dz)):
        for j in range(i + 1, len(dz)):
            gap = abs(dz[i] - dz[j])
            if min(gap, cols - gap) <= 1 and abs(dr[i] - dr[j]) <= 1:
                parent[root(i)] = root(j)
    group = np.full(len(cells), -1)
    group[np.nonzero(dense)[0]] = [root(i) for i in range(len(dz))]
    member = group[inv]
    prof = []
    for g in np.unique(member[member >= 0]):
        sel = member == g
        angle = 2 * math.pi * z0[sel] / model.lead            # (a circular mean: may straddle 0)
        zc = math.atan2(np.sin(angle).mean(), np.cos(angle).mean()) % (2 * math.pi) * model.lead / (2 * math.pi)
        prof.append((zc, float(np.median(rho[sel]))))
    if len(prof) < 2:
        return None
    prof = np.array(sorted(prof))
    prof = np.r_[prof, prof[:1] + [model.lead, 0]]          # closed over one lead
    trial = Helical(model.a, model.d, model.pitch, prof, True)
    off = np.abs(trial._profile_distance(z0, rho))
    if np.mean(off <= 0.02) < 0.8:
        return None                     # the corners don't follow one profile
    model, prof = _refine(model, prof, mesh.pts[corners[off <= 0.02]])
    rho, z0, _ = model.local(mesh.pts[corners])
    trial = Helical(model.a, model.d, model.pitch, prof, True)
    off = np.abs(trial._profile_distance(z0, rho))
    if np.mean(off <= 0.01) < 0.8:
        return None
    # which side is material: facets' normals point away from it
    X, N = mesh.fcent[facets], mesh.fn[facets]
    radial = X - model.a - np.outer((X - model.a) @ model.d, model.d)
    outward = np.einsum("ij,ij->i", radial, N) * mesh.farea[facets]
    inside = outward.sum() > 0                            # normals point away from the axis: a bolt
    return Helical(model.a, model.d, model.pitch, prof, inside), facets


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
