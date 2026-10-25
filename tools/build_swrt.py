"""Build the SWRT operator: keep K mirror-pair units of the table, one gain each. No WRT call is needed.

Rows and gains: on a set of practice spectra (make_synthetic.py), find K units and gains g so that sum_u g_u Psi_u(E) reproduces
the full-table sum sum_all Psi(E). The full-table sum is WRT (to single-precision round-off, check_table.py), so the gains carry
WRT's absolute units and the operator is used with c = 1. Each practice spectrum is scaled by 1 / |full sum| (Snl ~ E^3,
otherwise a few large seas dominate). Search: greedy from empty (each step tries every unit not yet kept, adds the one that
lowers the residual sum of squares most, refits all gains by ordinary least squares), then swap passes (replace one kept
unit at a time by the best other unit) until no swap lowers the residual. With a limit L on the offsets (di, dj) of the
kept rows (--offsets, or --cost and --line for the largest L at a given cost), the search is that of the paper: greedy under
the limit and swaps, then random sets within the limit and random replacements of three to six kept units, each followed by
swaps, in one run per search seed; the set with the smallest residual is kept. The residual is closed form from the Gram matrix,
which is summed over the seeds and cached per seed in the work directory, with a fingerprint of the table and the spectra:
a cache made with another table or other spectra is computed again.
Optional check, not applied (--c-check): c = <P, S_WRT> / <P, P> from one WRT call on one JONSWAP spectrum.
usage: python3 build_swrt.py --work DIR [--seeds 11,1,2,3,4,5,6,7,8,9] [--ks 10,20,30,50] [--c-check]
       python3 build_swrt.py --work DIR --ks 36 --offsets 19 [--starts 6 --kicks 30 --search-seeds 11,12,13]   (the paper's SWRT)
       python3 build_swrt.py --work DIR --ks 30,32,34,36,38 --cost C0,A,B,TDIA --line 4   (L at 4 times the cost of DIA)
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


def offsets_of(T, pair):
    """the offsets (di, dj) of k3 from k1 of the two rows of each unit (one offset if dj = 0, else two)"""
    return [frozenset((int(T["di"][r]), int(T["dj"][r])) for r in p) for p in pair]


def limits(Ks, offsets=None, cost=None, line=None):
    """the offset limit L of each K: the same L for every K (offsets), or the largest L at `line` times the cost of DIA in the
    cost model t = C0 + A K + B L (cost = C0, A, B, TDIA in one unit), at most 2 K; None for every K without a limit"""
    if offsets is None and cost is None: return {K: None for K in Ks}
    if cost is not None:
        c0, a, b, tdia = cost; L = {K: int(np.floor((line * tdia - c0 - a * K) / b + 1e-9)) for K in Ks}
    else: L = {K: int(offsets) for K in Ks}
    bad = [K for K in Ks if L[K] < 1]
    if bad: raise SystemExit(f"no offset is left for K = {bad} at this cost")
    return {K: min(L[K], 2 * K) for K in Ks}


def search(G, b, yy, K, L, LOC, rng, starts=6, kicks=30):
    """K units whose rows use at most L offsets (LOC: offsets_of), as in the paper: greedy under the limit and swaps, then
    `starts` random sets within the limit and `kicks` random replacements of three to six units of the best set, each followed
    by swaps. Returns the set with the smallest normalized residual e2 (sorted) and that e2."""
    G, b = G / yy, b / yy; N = len(b)
    loci = lambda S: len(frozenset().union(*[LOC[u] for u in S])) if S else 0
    e2 = lambda S: float(1 - b[S] @ np.linalg.solve(G[np.ix_(S, S)], b[S]))

    def greedy():
        S = []
        while len(S) < K:
            c = [q for q in range(N) if q not in S and loci(S + [q]) <= L]
            if not c: raise SystemExit(f"no set of K = {K} units uses at most {L} offsets")
            S.append(min(c, key=lambda q: e2(S + [q])))
        return S

    def swap(S):
        S = list(S)
        while True:
            best = e2(S)
            for pos in range(len(S)):
                rest = S[:pos] + S[pos + 1:]; c = [q for q in range(N) if q not in rest and loci(rest + [q]) <= L]
                q = min(c, key=lambda q: e2(rest + [q]))
                if e2(rest + [q]) < e2(S) - 1e-12: S = rest + [q]
            if e2(S) >= best - 1e-12: return S

    best = swap(greedy()); eb = e2(best)
    for _ in range(starts):                                          # random sets within the limit
        S = []
        for q in rng.permutation(N):
            if len(S) == K: break
            if loci(S + [int(q)]) <= L: S.append(int(q))
        if len(S) < K: continue
        S = swap(S)
        if e2(S) < eb - 1e-12: best, eb = S, e2(S)
    for _ in range(kicks):                                           # random replacements of three to six units
        S = list(best)
        for pos in rng.choice(K, min(int(rng.integers(3, 7)), K), replace=False):
            c = [q for q in range(N) if q not in S and loci(S[:pos] + S[pos + 1:] + [q]) <= L]
            if c: S[pos] = int(rng.choice(c))
        S = swap(S)
        if e2(S) < eb - 1e-12: best, eb = S, e2(S)
    return sorted(int(u) for u in best), float(eb)


def jonswap(freq, dirs, fp=0.14, hs=2.0, gamma=3.3, s2=20.0, th0=187.5):
    """JONSWAP in f, sigma_a/b = 0.07/0.09, times D(theta) ~ cos^s2((theta - th0)/2), scaled to Hs (the c check)"""
    f = freq; df = np.gradient(f); sg = np.where(f <= fp, 0.07, 0.09)
    S = f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gamma ** np.exp(-((f - fp) ** 2) / (2 * (sg * fp) ** 2))
    D = np.abs(np.cos(np.deg2rad(dirs - th0) / 2)) ** s2; D /= D.sum() * 2 * np.pi / len(dirs)
    E = S[:, None] * D[None]; return E * (hs ** 2 / 16) / (E.sum(1) * 2 * np.pi / len(dirs) * df).sum()


def build(work, seeds, Ks, c_check=False, json_name="build.json", L_of=None, starts=6, kicks=30, search_seeds=(11, 12, 13)):
    B = Builder(os.path.join(work, "table.npz")); T, NU, pair = B.T, B.NU, B.pair
    if max(Ks) > NU: raise SystemExit(f"K = {max(Ks)} is larger than the {NU} configurations of the table")
    G, b, yy, gram_s = 0, 0, 0.0, 0.0
    for sd in seeds:                                                 # pooled draws: sum of the per-draw Grams
        E = np.load(os.path.join(work, f"synthetic_s{sd}.npz"))["E"].astype(float)
        Gi, bi, yi, si = B.gram(E, os.path.join(work, f"gram_s{sd}.npz")); G, b, yy, gram_s = G + Gi, b + bi, yy + yi, gram_s + si
    L_of = L_of or {K: None for K in Ks}; free = tuple(K for K in Ks if L_of[K] is None); t0 = time.time()
    sel, path = select(G, b, yy, free) if free else ({}, [])
    LOC = offsets_of(T, pair); runs = {}
    for K in Ks:                                                     # the paper's search under an offset limit, best of the seeds
        if L_of[K] is None: continue
        res = [(sd,) + search(G, b, yy, K, L_of[K], LOC, np.random.default_rng(sd), starts, kicks) for sd in search_seeds]
        sd, keep, _ = min(res, key=lambda r: r[2]); keep = np.array(keep)
        g = np.linalg.solve(G[np.ix_(keep, keep)], b[keep]); sel[K] = (keep, g, 1 - (yy - b[keep] @ g) / yy, None)
        runs[K] = dict(offset_limit=L_of[K], offsets=len(frozenset().union(*[LOC[u] for u in keep])), search_seed=sd,
                       e2_by_search_seed={str(r[0]): r[2] for r in res}, starts=starts, kicks=kicks)
        print(f"K={K:2d}, at most {L_of[K]} offsets: e2 " + ", ".join(f"{r[2]:.6f} (seed {r[0]})" for r in res)
              + f"; kept seed {sd}, {runs[K]['offsets']} offsets", flush=True)
    sel = {K: sel[K] for K in Ks}; search_s = time.time() - t0
    out = {"seeds": list(seeds), "n_probe": 1500 * len(seeds), "probe_fit_r2_greedy_path": path, "units_total": NU,
           "rows_total": len(T["di"]), "gram_seconds": gram_s, "search_seconds": search_s}
    if c_check:                                                      # one WRT call (xnl_ref must be set up for this grid)
        import xnl_ref
        nd = swrt.n_dirs(T); Ej = jonswap(T["freq"], (360.0 / nd) * (np.arange(nd) + 0.5))
        Yj = xnl_ref.wrt(Ej).ravel(); Uj = B.unit_responses(Ej)[0]; out["c_check_jonswap"] = dict(fp=0.14, hs=2.0, gamma=3.3, s2=20)
    for K, (keep, g, r2, r2g) in sel.items():
        rows = pair[keep].ravel(); fn = os.path.join(work, f"swrt_K{K}.npz")
        c_chk = float(((g @ Uj[keep]) * Yj).sum() / ((g @ Uj[keep]) ** 2).sum()) if c_check else np.nan
        lim = runs.get(K, dict(offset_limit=-1, search_seed=-1))
        np.savez(fn, K=K, units=keep, rows=rows, g_unit=g, g_row=np.repeat(g, 2), c=1.0, c_check=c_chk, probes="synthetic",
                 seeds=np.array(seeds), probe_fit_r2=r2, na=swrt.n_dirs(T), nk=len(T["freq"]),
                 offset_limit=lim["offset_limit"], search_seed=lim["search_seed"],
                 **{k: T[k][rows] for k in swrt.GEOM + ("amem", "node")}, **{k: T[k] for k in swrt.CONST},
                 note="SWRT: rows of WRT's quadrature (2 per unit: row + mirror), gains fitted to the full-table sum on practice spectra, c = 1")
        out[f"K{K}"] = dict(file=os.path.basename(fn), units=keep.tolist(), gains=g.tolist(), c_check=c_chk if c_check else None, probe_fit_r2=r2,
                            greedy_only_r2=r2g, neg_gains=int((g < 0).sum()), configs=[(int(T["di"][pair[u, 0]]), int(abs(T["dj"][pair[u, 0]]))) for u in keep],
                            **({"search": runs[K]} if K in runs else {}))
        print(f"K={K:2d}: practice-fit R2 {r2:.4f}" + (f" (greedy alone {r2g:.4f})" if r2g is not None else f" (e2 {1 - r2:.6f})")
              + f", gains {g.min():.3g}..{g.max():.3g} (negative {(g < 0).sum()})"
              + (f", c check {c_chk:.4f}" if c_check else "") + f" -> {fn}", flush=True)
    json.dump(out, open(os.path.join(work, json_name), "w"), indent=1)
    print(f"search {search_s:.1f} s")
    return out


def add_search_args(ap):
    ap.add_argument("--offsets", type=int, default=None, help="limit on the offsets (di, dj) of the kept rows (the paper: 19 for K = 36)")
    ap.add_argument("--cost", default=None, help="C0,A,B,TDIA: cost model t = C0 + A K + B L of S_nl and the cost of DIA, for --line")
    ap.add_argument("--line", type=float, default=None, help="with --cost: L is the largest limit at this many times the cost of DIA")
    ap.add_argument("--starts", type=int, default=6, help="random starts of the search under a limit (default 6)")
    ap.add_argument("--kicks", type=int, default=30, help="random replacements of the search under a limit (default 30)")
    ap.add_argument("--search-seeds", default="11,12,13", help="one search per seed, the best is kept (default 11,12,13)")


def search_kw(a, Ks):
    if (a.cost is None) != (a.line is None): raise SystemExit("--cost and --line go together")
    if a.cost is not None and a.offsets is not None: raise SystemExit("give --offsets or --cost and --line, not both")
    cost = tuple(float(x) for x in a.cost.split(",")) if a.cost else None
    return dict(L_of=limits(Ks, a.offsets, cost, a.line), starts=a.starts, kicks=a.kicks,
                search_seeds=tuple(int(s) for s in a.search_seeds.split(",")))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", required=True); ap.add_argument("--seeds", default="11,1,2,3,4,5,6,7,8,9")
    ap.add_argument("--ks", default="10,20,30,50"); ap.add_argument("--c-check", action="store_true",
                                                                       help="one WRT call (needs the grid: --f1 --ratio --nf --ndir)")
    ap.add_argument("--f1", type=float); ap.add_argument("--ratio", type=float); ap.add_argument("--nf", type=int); ap.add_argument("--ndir", type=int)
    add_search_args(ap); a = ap.parse_args()
    if a.c_check and None in (a.f1, a.ratio, a.nf, a.ndir): ap.error("--c-check needs the grid: --f1 --ratio --nf --ndir")
    if a.c_check:
        import xnl_ref
        freq, _ = xnl_ref.grid(a.f1, a.ratio, a.nf, a.ndir); xnl_ref.init(freq, a.ndir, os.path.join(a.work, "xnl"))
    Ks = tuple(int(k) for k in a.ks.split(","))
    build(a.work, [int(s) for s in a.seeds.split(",")], Ks, a.c_check, **search_kw(a, Ks))
