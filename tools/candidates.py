"""Candidate sets of K configurations for a choice by model simulations (SWRT9 of the paper, its section 3d): the local optima
that the search under an offset limit reaches from many starts, and the N distinct sets with the smallest e2 among them.

With few configurations, sets with a similar e2 can differ widely in model simulations, so that e2 alone does not choose the set.
The starts are
  (a) for each limit L of --limits and each seed of --start-seeds, a random set within the limit (the units in the order of
      numpy's default_rng(seed).permutation, each added while the rows use at most L offsets), improved by swaps;
  (b) for each limit L of --search-limits and each seed of --search-seeds, the search of build_swrt.py (greedy, swaps, six
      random sets and 30 random replacements).
The defaults are the paper's: K = 9; (a) L = 7-12 with seeds 100-139; (b) L = 6-12 with seeds 11-20; 40 sets kept. The gains
are least squares on the practice spectra, as for every operator of build_swrt.py. Then tools/screen.py runs the candidates in
the six training simulations of single_point/ and retains the one with the fewest time steps.
Needs the work directory of make_operator.py for the grid (table.npz, synthetic_s*.npz and their Gram matrices).
usage: python3 tools/candidates.py --work work_40x24 [--k 9] [--jobs 8]
       ->  work_40x24/candidates_K9/: candidates.json, and c01.npz/.ww3op ... c40.npz/.ww3op in the order of e2"""
import os, sys, json, time, argparse
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import swrt, build_swrt, export_operator

S = {}   # state shared with the worker processes (fork)


def ints(s):
    """'7-12' or '11,12,13' or '6-8,10' -> list of ints"""
    out = []
    for part in s.split(","):
        a, _, b = part.partition("-"); out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


def offsets(units):
    return len(frozenset().union(*[S["LOC"][u] for u in units])) if len(units) else 0


def one(job):
    """one start: ('start', L, seed) a random set within the limit and swaps, or ('search', L, seed) the search of build_swrt.py"""
    kind, L, seed = job; G, b, yy, LOC, K = S["G"], S["b"], S["yy"], S["LOC"], S["K"]; rng = np.random.default_rng(seed)
    if kind == "search":
        units, e2 = build_swrt.search(G, b, yy, K, L, LOC, rng)
    else:
        start = []
        for q in rng.permutation(len(b)):
            if len(start) == K: break
            if offsets(start + [int(q)]) <= L: start.append(int(q))
        if len(start) < K: return None
        units, e2 = build_swrt.search(G, b, yy, K, L, LOC, rng, start=start)
    return dict(kind=kind, limit=L, seed=seed, units=units, e2=e2, offsets=offsets(units))


def save(fn, T, pair, keep, g, seeds, r2, limit, seed):
    """the operator file of build_swrt.py for one candidate"""
    rows = pair[keep].ravel()
    np.savez(fn, K=len(keep), units=keep, rows=rows, g_unit=g, g_row=np.repeat(g, 2), c=1.0, c_check=np.nan, probes="synthetic",
             seeds=np.array(seeds), probe_fit_r2=r2, na=swrt.n_dirs(T), nk=len(T["freq"]), offset_limit=limit, search_seed=seed,
             **{k: T[k][rows] for k in swrt.GEOM + ("amem", "node")}, **{k: T[k] for k in swrt.CONST},
             note="SWRT candidate of candidates.py: rows of WRT's quadrature (2 per unit: row + mirror), gains fitted to the "
                  "full-table sum on practice spectra, c = 1")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", required=True, help="work directory of make_operator.py for the grid")
    ap.add_argument("--k", type=int, default=9, help="number of configurations (default 9, SWRT9)")
    ap.add_argument("--limits", default="7-12", help="offset limits of the random starts (a) (default 7-12)")
    ap.add_argument("--start-seeds", default="100-139", help="seeds of the random starts (a) (default 100-139)")
    ap.add_argument("--search-limits", default="6-12", help="offset limits of the searches (b) (default 6-12)")
    ap.add_argument("--search-seeds", default="11-20", help="seeds of the searches (b) (default 11-20; 'none' for no search)")
    ap.add_argument("--keep", type=int, default=40, help="number of distinct sets kept, smallest e2 first (default 40)")
    ap.add_argument("--seeds", default="11,1,2,3,4,5,6,7,8,9", help="practice-spectrum seeds (default: the paper's ten)")
    ap.add_argument("--jobs", type=int, default=os.cpu_count(), help="processes (default: all cores)")
    a = ap.parse_args(); t0 = time.time(); work = os.path.abspath(a.work); K = a.k; seeds = ints(a.seeds)
    B = build_swrt.Builder(os.path.join(work, "table.npz")); G, b, yy = 0, 0, 0.0
    for sd in seeds:                                                 # the Gram matrices of make_operator.py (computed if missing)
        E = np.load(os.path.join(work, f"synthetic_s{sd}.npz"))["E"].astype(float)
        Gi, bi, yi, _ = B.gram(E, os.path.join(work, f"gram_s{sd}.npz")); G, b, yy = G + Gi, b + bi, yy + yi
    S.update(G=G, b=b, yy=yy, LOC=build_swrt.offsets_of(B.T, B.pair), K=K)
    jobs = [("start", L, s) for L in ints(a.limits) for s in ints(a.start_seeds)]
    if a.search_seeds != "none": jobs += [("search", L, s) for L in ints(a.search_limits) for s in ints(a.search_seeds)]
    with Pool(a.jobs) as p: res = [r for r in p.map(one, jobs) if r]
    best = {}
    for r in sorted(res, key=lambda r: r["e2"]): best.setdefault(tuple(r["units"]), r)   # distinct sets, smallest e2 first
    keep = list(best.values())[:a.keep]
    out = os.path.join(work, f"candidates_K{K}"); os.makedirs(out, exist_ok=True); cands = []
    for i, r in enumerate(keep, 1):
        u = np.array(r["units"]); g = np.linalg.solve(G[np.ix_(u, u)], b[u]); r2 = float(1 - (yy - b[u] @ g) / yy)
        name = f"c{i:02d}"; npz = os.path.join(out, f"{name}.npz"); save(npz, B.T, B.pair, u, g, seeds, r2, r["limit"], r["seed"])
        cands.append(dict(name=name, units=r["units"], e2=r["e2"], offsets=r["offsets"], start=r["kind"], limit=r["limit"],
                          seed=r["seed"], gains=g.tolist(), npz=os.path.basename(npz),
                          ww3op=os.path.basename(export_operator.export(npz))))
    json.dump(dict(K=K, runs=len(jobs), local_optima=len(res), distinct=len(best), kept=len(keep), seeds=seeds,
                   args=dict(limits=a.limits, start_seeds=a.start_seeds, search_limits=a.search_limits, search_seeds=a.search_seeds),
                   candidates=cands, all=res, minutes=(time.time() - t0) / 60), open(os.path.join(out, "candidates.json"), "w"), indent=1)
    print(f"{len(jobs)} starts, {len(best)} distinct local optima; kept {len(keep)} with e2 {keep[0]['e2']:.4f} to "
          f"{keep[-1]['e2']:.4f} -> {out}/candidates.json ({(time.time() - t0) / 60:.0f} min)")


if __name__ == "__main__":
    main()
