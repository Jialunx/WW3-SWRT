"""SWRT at run time: WRT's own pair sum, restricted to a few rows of its quadrature, one gain per mirror-pair unit, on any grid
of the table (nk frequencies with a constant ratio, na directions).

For one row (di, dj, u2, v2, u4, v4, zz(i)) and every anchor cell (i, j), exactly as Xnl's q_xnl4v4 / q_t13v4:
    N    = E / (2 pi sigma) * cg / k                                   action density on the (kx, ky) plane
    N1   = N(i, j),   N3 = N(i + di, j + dj),   N2, N4 bilinear at (i + u2, j + v2), (i + u4, j + v4)
           (frequency corners above the highest bin are clamped to it; directions wrap)
    t    = zz(i) * [N1 N3 (N4 - N2) + N2 N4 (N3 - N1)]                  skipped if N1 and N3 < qf_frac Nmax (k1/k3)^7.5
    dS(i, j)            += 2 t * a(i3) * b(i)                          a = k dtheta dk   (area of the partner cell)
    dS(i + di, j + dj)  -= 2 t * a(i)  * b(i3)                         b = k / cg * 2 pi sigma  (N -> E(f, theta))
The anchor gains what the partner loses in action (a(i3) a(i) on both sides), so every row conserves wave action.
SWRT = c * sum over kept units g_u * (row + its mirror row). The full table (all units, g = 1, c = 1) is WRT.

    from swrt import SWRT
    S = SWRT("swrt_K30.npz").snl(E)      # E (nk, na) or (n, nk, na)
"""
import numpy as np

TPI = 2 * np.pi
GEOM = ("di", "dj", "u2", "v2", "u4", "v4", "zz")
CONST = ("freq", "sigma", "k", "cg", "a", "b", "scale", "qf_frac", "kpow")


def load_table(path):
    return dict(np.load(path))


def n_dirs(T):
    """number of directions of a table or operator (files of the paper's scripts do not store it: 24)"""
    return int(T["na"]) if "na" in T else 24


class Rows:
    """gather/scatter plan for a set of rows; responses() gives one map per row, snl() the weighted sum"""

    def __init__(self, T, rows=None):
        rows = np.arange(len(T["di"])) if rows is None else np.asarray(rows)
        self.c = {k: np.asarray(T[k], float) for k in CONST}
        self.nk, self.na = NK, NA = len(self.c["freq"]), n_dirs(T)
        self.di, self.dj = T["di"][rows].astype(int), T["dj"][rows].astype(int)
        self.zz = np.asarray(T["zz"][rows], float); R = len(rows); I = np.arange(NK); J = np.arange(NA)
        i3 = I[None] + self.di[:, None]; self.valid = (i3 <= NK - 1).astype(float); self.i3 = np.minimum(i3, NK - 1)   # (R, nk)
        self.idx3 = (self.i3[:, :, None] * NA + (J[None, None] + self.dj[:, None, None]) % NA)                        # (R, nk, na)
        self.corners = [self._corners(np.asarray(T[u][rows], float), np.asarray(T[v][rows], float)) for u, v in (("u2", "v2"), ("u4", "v4"))]
        self.fa = self.c["a"][self.i3] * self.c["b"][None, :]            # anchor factor a(i3) b(i)   (R, nk)
        self.fp = self.c["a"][None, :] * self.c["b"][self.i3]            # partner factor a(i) b(i3)
        self.qpow = self.c["kpow"][self.i3]                              # (k1 / k3)^7.5 density threshold shape
        self.cfg = {}
        for r, key in enumerate(zip(self.di, self.dj)): self.cfg.setdefault(key, []).append(r)

    def _corners(self, U, V):
        NK, NA = self.nk, self.na
        I = np.arange(NK)[None]; J = np.arange(NA)[None, None]
        i0 = I + np.floor(U)[:, None].astype(int); fu = np.broadcast_to((U - np.floor(U))[:, None], i0.shape).copy()
        lo, hi = i0 < 0, i0 > NK - 1; fu[lo | hi] = 0.0; i0 = np.clip(i0, 0, NK - 1); ip = np.minimum(i0 + 1, NK - 1)
        j0 = np.floor(V).astype(int)[:, None, None]; fv = (V - np.floor(V))[:, None, None]
        ja, jb = (J + j0) % NA, (J + j0 + 1) % NA; fu = fu[:, :, None]
        idx = [i0[:, :, None] * NA + ja, ip[:, :, None] * NA + ja, i0[:, :, None] * NA + jb, ip[:, :, None] * NA + jb]
        w = [(1 - fu) * (1 - fv), fu * (1 - fv), (1 - fu) * fv, fu * fv]
        return idx, w

    def _n(self, E):
        """densities of every row at every anchor cell, n1, n2, n3, n4 (B, R, nk, na) (n1 broadcast), and the density filter"""
        E = np.asarray(E, float); E = E[None] if E.ndim == 2 else E; B = len(E)
        N = (E / (TPI * self.c["sigma"])[None, :, None] * (self.c["cg"] / self.c["k"])[None, :, None]).reshape(B, -1)
        n1 = N.reshape(B, 1, self.nk, self.na); n3 = N[:, self.idx3]
        n2 = sum(w[None] * N[:, i] for i, w in zip(*self.corners[0])); n4 = sum(w[None] * N[:, i] for i, w in zip(*self.corners[1]))
        qmin = (self.c["qf_frac"] * N.max(1))[:, None, None, None] * self.qpow[None, :, :, None]
        return n1, n2, n3, n4, (n1 < qmin) & (n3 < qmin), B

    def _t(self, E):
        """per-row transfer t (B, R, nk, na) in N units, after the density filter"""
        n1, n2, n3, n4, skip, B = self._n(E)
        t = (self.zz * self.valid)[None, :, :, None] * (n1 * n3 * (n4 - n2) + n2 * n4 * (n3 - n1))
        return np.where(skip, 0.0, t) * self.c["scale"], B

    def _scatter(self, ta, tp):
        """ta (B, R, nk, na) deposited with the anchor factor at (i, j), minus tp with the partner factor at (i + di, j + dj)"""
        out = np.einsum("brij,ri->bij", ta, self.fa)
        for (di, dj), rs in self.cfg.items():                  # partner deposits, one shifted add per (di, dj)
            part = np.einsum("brij,ri->bij", tp[:, rs], self.fp[rs])
            out[:, di:, :] -= np.roll(part[:, :self.nk - di, :], dj, axis=-1)
        return out

    def responses(self, E):
        """(B, R, nk, na): each row's own Snl map, with WRT's factors, in E(f, theta) units"""
        t, B = self._t(E); out = t * self.fa[None, :, :, None]; part = t * self.fp[None, :, :, None]
        for (di, dj), rs in self.cfg.items():
            out[:, rs, di:, :] -= np.roll(part[:, rs, :self.nk - di, :], dj, axis=-1)
        return out

    def snl(self, E, g=None):
        """sum over rows with per-row gains g (default 1)"""
        t, B = self._t(E); g = np.ones(len(self.di)) if g is None else np.asarray(g, float)
        t = t * g[None, :, None, None]
        return self._scatter(t, t)

    def diag(self, E, g=None):
        """diagonal term dS(f, theta) / dE(f, theta) of the same sum, as Xnl computes it for WRT (q_t13v4, diagk1 / diagk3):
        every term differentiated with respect to its own N1 (deposited at the anchor) and N3 (at the partner), with N2 and N4,
        interpolated on the locus, held fixed. With the full table and g = 1 this is Xnl's diagonal."""
        n1, n2, n3, n4, skip, B = self._n(E); g = np.ones(len(self.di)) if g is None else np.asarray(g, float)
        z = (self.zz * self.valid * g[:, None])[None, :, :, None] * self.c["scale"]
        d1 = np.where(skip, 0.0, z * (n3 * (n4 - n2) - n2 * n4)); d3 = np.where(skip, 0.0, z * (n1 * (n4 - n2) + n2 * n4))
        return self._scatter(d1, d3) * (self.c["cg"] / (TPI * self.c["sigma"] * self.c["k"]))[None, :, None]   # x dN/dE


class SWRT:
    """deploy object for an operator file written by build_swrt.py (rows, per-row gains, c, run-time constants)"""

    def __init__(self, path, use_c=True):
        d = np.load(path); self.K = int(d["K"]); self.c = float(d["c"]) if use_c else 1.0
        self.rows = Rows({k: d[k] for k in GEOM + CONST + (("na",) if "na" in d else ())}); self.g = d["g_row"].astype(float)

    def snl(self, E):
        E = np.asarray(E, float); one = E.ndim == 2
        S = self.c * self.rows.snl(E, self.g)
        return S[0] if one else S

    def diag(self, E):
        """diagonal term dS / dE per bin (Rows.diag), as WAVEWATCH III's NL2 (Xnl) returns it for implicit time stepping"""
        E = np.asarray(E, float); one = E.ndim == 2
        D = self.c * self.rows.diag(E, self.g)
        return D[0] if one else D
