"""Decode WRT's quadruplet database quad99999.bqf (Xnl 5.03, Fortran sequential unformatted) of any grid, in either byte
order (WAVEWATCH III writes it big-endian by default, the tools' Xnl build in the machine's order).

Record layout (written by q_makegrid in mod_xnl4v5.f90; nk frequencies, na directions, nl points per locus):
   1  q_header (21 chars, '(3i3.3,6i2.2)': naq nkq nlocus0 iq_grid iq_geom iq_disp iq_cple iq_locus iq_interp)
   2  naq, nkq                      3  q_sig (nk)              4  q_ad (na, degrees)
   5  iq_geom, iq_disp, iq_geom     6  q_depth                 7  quad_nloc (nk (nk+1) / 2, na / 2 + 1)
   8-11  quad_ik2, ia2, ik4, ia4 (nk (nk+1) / 2, na / 2 + 1, nl)   12-19  quad_w1k2..w4k2, w1k4..w4k4   20  quad_zz
Indexing: kmem = jk3 + (jk1-1)*nk - jk1(jk1-1)/2 (jk3 >= jk1, 1-based); amem = |ia3 - ia1| + 1 (0..180 deg);
directions relative to ia1 = 1, may be negative (wrap mod na). Bilinear corners w1 (ik,ia), w2 (ik+1,ia),
w3 (ik,ia+1), w4 (ik+1,ia+1); the frequency fraction is linear in k (q_weight). ik <= nk; a corner ik+1 = nk+1 is
clamped to nk at run time and the f^-5 tail is already inside zz (tail factor wt_k2*wt_k4).
Slots il > nloc are never read by WRT and hold uninitialised values (they differ between regenerated copies)."""
import numpy as np

NAMES = ["ik2", "ia2", "ik4", "ia4", "w1k2", "w2k2", "w3k2", "w4k2", "w1k4", "w2k4", "w3k4", "w4k4", "zz"]


def decode(path):
    raw = open(path, "rb").read(); recs = []; p = 0
    bo = "<" if int(np.frombuffer(raw[:4], "<i4")[0]) == 21 else ">"         # the header record holds 21 characters
    while p < len(raw):
        n = int(np.frombuffer(raw[p:p + 4], bo + "i4")[0]); recs.append((p + 4, n)); p += 4 + n + 4
    arr = lambda k, dt: np.frombuffer(raw[recs[k][0]:recs[k][0] + recs[k][1]], bo + dt[1:])
    na, nk = (int(x) for x in arr(1, "<i4")[:2])
    mk, ma = nk * (nk + 1) // 2, na // 2 + 1
    if len(recs) != 20 or recs[6][1] != 4 * mk * ma or recs[7][1] % (4 * mk * ma):
        raise ValueError(f"{path}: not a full-storage database of {nk} frequencies and {na} directions (records {len(recs)})")
    nl = recs[7][1] // (4 * mk * ma)
    out = {"header": raw[recs[0][0]:recs[0][0] + recs[0][1]].decode(), "byte_order": bo, "naq_nkq": arr(1, "<i4").copy(), "nk": nk, "na": na, "nl": nl,
           "sigma": arr(2, "<f4").astype(float), "angle": arr(3, "<f4").astype(float), "flags": arr(4, "<i4").copy(),
           "depth": float(arr(5, "<f4")[0]), "nloc": arr(6, "<i4").reshape((mk, ma), order="F"),
           "nrec": len(recs), "rec_bytes": [n for _, n in recs], "file_bytes": len(raw)}
    for i, nm in enumerate(NAMES):
        out[nm] = arr(7 + i, "<i4" if i < 4 else "<f4").reshape((mk, ma, nl), order="F")
    return out


def kmem(j1, j3, nk):
    """0-based slot of the 1-based pair (jk1, jk3 >= jk1)"""
    return j3 + (j1 - 1) * nk - j1 * (j1 - 1) // 2 - 1


def raw_nodes(Q, j1, j3, am):
    """the 13 stored numbers of every node of pair (jk1, jk3 >= jk1), amem am (1-based), as an (n, 13) array"""
    m = kmem(j1, j3, Q["nk"]); n = int(Q["nloc"][m, am - 1])
    return np.column_stack([np.asarray(Q[k][m, am - 1, :n], float) for k in NAMES]) if n > 0 else np.zeros((0, 13))


def offsets(Q, j1, j3, am):
    """(u2, v2, u4, v4, zz) per node: k2/k4 positions relative to (jk1, ia1) in bins, fractional parts from the weights"""
    r = raw_nodes(Q, j1, j3, am)
    if not len(r): return np.zeros((0, 5))
    u2 = r[:, 0] + r[:, 5] + r[:, 7] - j1; v2 = r[:, 1] - 1 + r[:, 6] + r[:, 7]
    u4 = r[:, 2] + r[:, 9] + r[:, 11] - j1; v4 = r[:, 3] - 1 + r[:, 10] + r[:, 11]
    return np.column_stack([u2, v2, u4, v4, r[:, 12]])
