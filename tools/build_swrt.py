"""Build the SWRT operator: keep K mirror-pair units of the table, one gain each. No WRT call is needed.

Rows and gains: on a set of practice spectra (make_synthetic.py), find K units and gains g so that sum_u g_u Psi_u(E) reproduces
the full-table sum sum_all Psi(E). The full-table sum is WRT (to single-precision round-off, check_table.py), so the gains carry
WRT's absolute units and the operator is used with c = 1. Each practice spectrum is scaled by 1 / |full sum| (Snl ~ E^3,
otherwise a few large seas dominate). Search: greedy from empty (each step tries every unit not yet kept, adds the one that
lowers the residual sum of squares most, refits all gains by ordinary least squares), then swap passes (replace one kept
unit at a time by the best other unit) until no swap lowers the residual. The residual is closed form from the Gram matrix,
which is summed over the seeds and cached per seed in the work directory, with a fingerprint of the table and the spectra:
a cache made with another table or other spectra is computed again.
Optional check, not applied (--c-check): c = <P, S_WRT> / <P, P> from one WRT call on one JONSWAP spectrum.
usage: python3 build_swrt.py --work DIR [--seeds 11,1,2,3,4,5,6,7,8,9] [--ks 10,20,30,50] [--c-check]
       ->  DIR/swrt_K{K}.npz, DIR/build.json   (needs DIR/table.npz and DIR/synthetic_s{seed}.npz)"""
import os, json, time, hashlib, argparse
import numpy as np
import swrt

BATCH = 16
FP_KEYS = swrt.GEOM + swrt.CONST + ("unit",)


def fingerprint(T, E):
    """sha256 of the table arrays that enter the responses and of the practice spectra"""
    h = hashlib.sha256()
    for k in FP_KEYS: h.update(k.encode()); h.update(np.ascontiguousarray(T[k], float).tobytes())
    h.update(str(swrt.n_dirs(T)).encode()); h.update(np.ascontiguousarray(E, float).tobytes())
    return h.hexdigest()


class Builder:
    def __init__(self, table):
        self.T = swrt.load_table(table); self.R = swrt.Rows(self.T); self.NU = NU = int(self.T["unit"].max()) + 1
        self.pair = np.array([np.flatnonzero(self.T["unit"] == u) for u in range(NU)])       # (NU, 2) row ids of each unit

    def unit_responses(self, E):
        r = self.R.responses(E); r = r.reshape(r.shape[0], r.shape[1], -1)
        return r[:, self.pair[:, 0]] + r[:, self.pair[:, 1]]                                 # (B, NU, nk na)

    def gram(self, E, cache):
        """normalised Gram matrix of the unit responses and their projections on the full-table sum (cached)"""
        fp = fingerprint(self.T, E)
        if os.path.exists(cache):
            z = np.load(cache)
            if "fp" in z and str(z["fp"]) == fp: return z["G"], z["b"], float(z["yy"]), float(z["seconds"])
            print(f"{os.path.basename(cache)} was made with another table or other spectra: computing it again", flush=True)
        NU = self.NU; G = np.zeros((NU, NU)); b = np.zeros(NU); yy = 0.0; t0 = time.time()
        for s in range(0, len(E), BATCH):
            U = self.unit_responses(E[s:s + BATCH]); full = U.sum(1); sc = 1 / np.sqrt((full ** 2).sum(1))
            A = (U * sc[:, None, None]).transpose(1, 0, 2).reshape(NU, -1); y = (full * sc[:, None]).ravel()
            G += A @ A.T; b += A @ y; yy += float(y @ y)
        sec = time.time() - t0; np.savez(cache, G=G, b=b, yy=yy, seconds=sec, fp=fp)
        print(f"Gram over {len(E)} practice spectra, {NU} units: {sec:.0f} s", flush=True)
        return G, b, yy, sec


def select(G, b, yy, Ks):
    NU = len(b)
    fit = lambda Ss: np.linalg.solve(G[np.ix_(Ss, Ss)], b[Ss]); sse = lambda Ss: yy - b[Ss] @ fit(Ss)   # OLS residual
    S, path = [], []
    while len(S) < max(Ks):                                         # greedy over all units
        S.append(min((q for q in range(NU) if q not in S), key=lambda q: sse(S + [q]))); path.append(1 - sse(S) / yy)

    def swap(S):
        S = list(S)
        while True:
            best = sse(S)
            for pos in range(len(S)):
                rest = S[:pos] + S[pos + 1:]; q = min((q for q in range(NU) if q not in rest), key=lambda q: sse(rest + [q]))
                if sse(rest + [q]) < sse(S) - 1e-12 * yy: S = rest + [q]
            if sse(S) >= best - 1e-12 * yy: return S
    out = {}
    for K in Ks:
        keep = swap(S[:K]); out[K] = (np.array(keep), fit(keep), 1 - sse(keep) / yy, path[K - 1])
    return out, path


def jonswap(freq, dirs, fp=0.14, hs=2.0, gamma=3.3, s2=20.0, th0=187.5):
    """JONSWAP in f, sigma_a/b = 0.07/0.09, times D(theta) ~ cos^s2((theta - th0)/2), scaled to Hs (the c check)"""
    f = freq; df = np.gradient(f); sg = np.where(f <= fp, 0.07, 0.09)
    S = f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gamma ** np.exp(-((f - fp) ** 2) / (2 * (sg * fp) ** 2))
    D = np.abs(np.cos(np.deg2rad(dirs - th0) / 2)) ** s2; D /= D.sum() * 2 * np.pi / len(dirs)
    E = S[:, None] * D[None]; return E * (hs ** 2 / 16) / (E.sum(1) * 2 * np.pi / len(dirs) * df).sum()


def build(work, seeds, Ks, c_check=False, json_name="build.json"):
    B = Builder(os.path.join(work, "table.npz")); T, NU, pair = B.T, B.NU, B.pair
    if max(Ks) > NU: raise SystemExit(f"K = {max(Ks)} is larger than the {NU} configurations of the table")
    G, b, yy, gram_s = 0, 0, 0.0, 0.0
    for sd in seeds:                                                 # pooled draws: sum of the per-draw Grams
        E = np.load(os.path.join(work, f"synthetic_s{sd}.npz"))["E"].astype(float)
        Gi, bi, yi, si = B.gram(E, os.path.join(work, f"gram_s{sd}.npz")); G, b, yy, gram_s = G + Gi, b + bi, yy + yi, gram_s + si
    t0 = time.time(); sel, path = select(G, b, yy, Ks); search_s = time.time() - t0
    out = {"seeds": list(seeds), "n_probe": 1500 * len(seeds), "probe_fit_r2_greedy_path": path, "units_total": NU,
           "rows_total": len(T["di"]), "gram_seconds": gram_s, "search_seconds": search_s}
    if c_check:                                                      # one WRT call (xnl_ref must be set up for this grid)
        import xnl_ref
        nd = swrt.n_dirs(T); Ej = jonswap(T["freq"], (360.0 / nd) * (np.arange(nd) + 0.5))
        Yj = xnl_ref.wrt(Ej).ravel(); Uj = B.unit_responses(Ej)[0]; out["c_check_jonswap"] = dict(fp=0.14, hs=2.0, gamma=3.3, s2=20)
    for K, (keep, g, r2, r2g) in sel.items():
        rows = pair[keep].ravel(); fn = os.path.join(work, f"swrt_K{K}.npz")
        c_chk = float(((g @ Uj[keep]) * Yj).sum() / ((g @ Uj[keep]) ** 2).sum()) if c_check else np.nan
        np.savez(fn, K=K, units=keep, rows=rows, g_unit=g, g_row=np.repeat(g, 2), c=1.0, c_check=c_chk, probes="synthetic",
                 seeds=np.array(seeds), probe_fit_r2=r2, na=swrt.n_dirs(T), nk=len(T["freq"]),
                 **{k: T[k][rows] for k in swrt.GEOM + ("amem", "node")}, **{k: T[k] for k in swrt.CONST},
                 note="SWRT: rows of WRT's quadrature (2 per unit: row + mirror), gains fitted to the full-table sum on practice spectra, c = 1")
        out[f"K{K}"] = dict(file=os.path.basename(fn), units=keep.tolist(), gains=g.tolist(), c_check=c_chk, probe_fit_r2=r2,
                            greedy_only_r2=r2g, neg_gains=int((g < 0).sum()), configs=[(int(T["di"][pair[u, 0]]), int(abs(T["dj"][pair[u, 0]]))) for u in keep])
        print(f"K={K:2d}: practice-fit R2 {r2:.4f} (greedy alone {r2g:.4f}), gains {g.min():.3g}..{g.max():.3g} (negative {(g < 0).sum()})"
              + (f", c check {c_chk:.4f}" if c_check else "") + f" -> {fn}", flush=True)
    json.dump(out, open(os.path.join(work, json_name), "w"), indent=1)
    print(f"search {search_s:.1f} s")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", required=True); ap.add_argument("--seeds", default="11,1,2,3,4,5,6,7,8,9")
    ap.add_argument("--ks", default="10,20,30,50"); ap.add_argument("--c-check", action="store_true",
                                                                       help="one WRT call (needs the grid: --f1 --ratio --nf --ndir)")
    ap.add_argument("--f1", type=float); ap.add_argument("--ratio", type=float); ap.add_argument("--nf", type=int); ap.add_argument("--ndir", type=int)
    a = ap.parse_args()
    if a.c_check and None in (a.f1, a.ratio, a.nf, a.ndir): ap.error("--c-check needs the grid: --f1 --ratio --nf --ndir")
    if a.c_check:
        import xnl_ref
        freq, _ = xnl_ref.grid(a.f1, a.ratio, a.nf, a.ndir); xnl_ref.init(freq, a.ndir, os.path.join(a.work, "xnl"))
    build(a.work, [int(s) for s in a.seeds.split(",")], tuple(int(k) for k in a.ks.split(",")), a.c_check)
