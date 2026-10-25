"""Write an SWRT operator (swrt_K*.npz of build_swrt.py) as the plain-text file that WAVEWATCH III's NL6 reads
(ww3_src/model/src/w3snl6md.F90). The file holds the row arrays of the compiled kernel, so that W3SNL6 only reads and applies them.

usage: python3 export_operator.py --op DIR/swrt_K36.npz [--out FILE]   ->  FILE (default DIR/swrt_K36.ww3op)

Format (version 1, plain text, read with list-directed Fortran READ; integers %d, reals %.9e, which round-trips REAL*4):
    WW3-SWRT-OPERATOR          magic line
    1                          format version
    nk na nr K                 frequencies, directions, rows, mirror-pair units (nr = 2 K)
    swrt_K36.npz               source operator file (name only)
followed by one block per array: a line with the array name, then its values (fastest index first):
    sigma(nk)                  radian frequencies of the table grid (rad/s); W3SNL6 checks WAVEWATCH III's SIG against them
    di dj(nr)                  k3 offset in frequency / direction bins
    i2 f2 j2 g2(nr)            floor and fraction of the k2 offset: frequency (i2, f2) and direction (j2, g2); i4 f4 j4 g4 for k4
    w fa fp(nk, nr)            zz * gain * scale (0 where k3 is off the grid); a(i3) b(i1); a(i1) b(i3)
    conv(nk)                   E(f, theta) -> N(kx, ky): cg / (2 pi sigma k)
    qfrac                      density filter: skip a pair if N1 and N3 < qfrac * max(N) * kpow(i3)
    kpow(nk)                   (k1 / k)^7.5
The gains are the operator's g_row times its c (c = 1 for every operator written by build_swrt.py).
The rows are written grouped by resonance curve (sorted by di, then dj; stable), because the default kernel of W3SNL6 shares the
k1 and k3 work of the rows of one curve. The order of the rows changes the result only by round-off.
W3SNL6 stops with a message unless nk = NK, na = NTH and sigma = SIG to a relative 1e-5."""
import os, argparse
import numpy as np
import swrt

MAGIC, VERSION, TPI = "WW3-SWRT-OPERATOR", 1, 2 * np.pi
NAMES = ("di", "dj", "i2", "f2", "j2", "g2", "i4", "f4", "j4", "g4", "w", "fa", "fp")


def rows_of(d, g_row, T):
    """the kernel's row arrays from the operator's rows and gains and its run-time constants T"""
    NK = len(T["sigma"])
    R = len(d["di"]); di = np.asarray(d["di"], int); I = np.arange(NK)
    i3 = I[None] + di[:, None]; ok = i3 <= NK - 1; i3c = np.minimum(i3, NK - 1)
    w = np.where(ok, np.asarray(d["zz"], float) * np.asarray(g_row, float)[:, None] * float(T["scale"]), 0.0)
    fa, fp = T["a"][i3c] * T["b"][None, :], T["a"][None, :] * T["b"][i3c]
    fl = lambda x: np.floor(np.asarray(x, float)).astype(np.int32); fr = lambda x: (np.asarray(x, float) - np.floor(np.asarray(x, float))).astype(np.float32)
    F = lambda x: np.asfortranarray(x.T, dtype=np.float32)                     # (nk, R)
    return (di.astype(np.int32), np.asarray(d["dj"], np.int32), fl(d["u2"]), fr(d["u2"]), fl(d["v2"]), fr(d["v2"]),
            fl(d["u4"]), fr(d["u4"]), fl(d["v4"]), fr(d["v4"]), F(w), F(fa), F(fp))


def kernel_arrays(op):
    """all arrays of the file, from the operator's rows, gains and constants"""
    d = np.load(op); T = {k: np.asarray(d[k]) for k in swrt.CONST}
    sig, k, cg = T["sigma"], T["k"], T["cg"]
    out = dict(zip(NAMES, rows_of(d, np.asarray(d["g_row"], float) * float(d["c"]), T)))
    out.update(sigma=np.asarray(sig, float), conv=(cg / (TPI * sig * k)).astype(np.float32), qfrac=np.float32(T["qf_frac"]),
               kpow=T["kpow"].astype(np.float32), K=int(d["K"]), nk=len(sig), na=swrt.n_dirs(d), nr=len(d["di"]), source=os.path.basename(op))
    return out


def write(path, a, per_line=8):
    def block(f, name, v, fmt):
        f.write(name + "\n"); v = np.ravel(v, order="F")
        for s in range(0, len(v), per_line): f.write(" ".join(fmt % x for x in v[s:s + per_line]) + "\n")
    with open(path, "w") as f:
        f.write(f"{MAGIC}\n{VERSION}\n{a['nk']} {a['na']} {a['nr']} {a['K']}\n{a['source']}\n")
        block(f, "sigma", a["sigma"], "%.9e")
        for n in NAMES:
            v = a[n]; assert np.all(np.isfinite(v)), n
            block(f, n, v, "%d" if v.dtype.kind == "i" else "%.9e")
        block(f, "conv", a["conv"], "%.9e"); block(f, "qfrac", [a["qfrac"]], "%.9e"); block(f, "kpow", a["kpow"], "%.9e")


def group_rows(a):
    """the rows of one resonance curve (same di, dj) next to each other, sorted by di, then dj (stable)"""
    o = np.lexsort((a["dj"], a["di"]))
    for n in NAMES: a[n] = a[n][o] if a[n].ndim == 1 else np.asfortranarray(a[n][:, o])
    return a


def export(op, out=None):
    out = out or os.path.splitext(op)[0] + ".ww3op"
    arr = group_rows(kernel_arrays(op)); write(out, arr)
    print(f"wrote {out}: K = {arr['K']}, {arr['nr']} rows, nk = {arr['nk']}, na = {arr['na']} (from {arr['source']})")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--op", required=True, help="operator .npz of build_swrt.py"); ap.add_argument("--out", default=None)
    a = ap.parse_args(); export(a.op, a.out)
