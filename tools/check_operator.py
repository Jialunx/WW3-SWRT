"""Checks against WRT on held-out practice spectra of the grid (make_synthetic.py with a seed not used for the selection,
default 1000; every --stride-th of its 1500 spectra). The reference is the deep-water WRT of xnl_ref.wrt().
  table      the full table summed with WRT's run-time rules (swrt.Rows) against WRT: max |P - Y| / max |Y| per spectrum.
             On the paper's grid this is below 5e-5 (single-precision round-off).
  operators  each operator's error, eps of van Vledder (2006, Eq. 58) on the direction-integrated transfer:
             sum_f |P(f) - Y(f)| df / sum_f |Y(f)| df, df = f (r - 1/r) / 2; median by spectrum type and overall.
The practice-spectrum ranges are those of the work directory (grid.json of make_operator.py), or the defaults.
usage: python3 check_operator.py --f1 0.03453 --ratio 1.1 --nf 40 --ndir 24 --work DIR [--ops DIR/swrt_K36.npz ...]
       ->  DIR/check.json   (needs DIR/table.npz; makes DIR/synthetic_s1000.npz if missing)"""
import os, json, argparse
import numpy as np
import swrt, xnl_ref, make_synthetic

TABLE_TOL = 1e-3          # a table that differs more from WRT is wrong


def family(work):
    T = swrt.load_table(os.path.join(work, "table.npz")); g = os.path.join(work, "grid.json")
    if os.path.exists(g):
        d = json.load(open(g)); return make_synthetic.family(T, tuple(d["fp_range"]), tuple(d["swell_range"]))
    return make_synthetic.family(T)


def heldout(work, seed=1000, stride=5, fam=None):
    """held-out practice spectra (made, or made again, for this grid and these ranges), their labels and WRT's transfer"""
    fam = fam or family(work); fn = os.path.join(work, f"synthetic_s{seed}.npz")
    if not (os.path.exists(fn) and fam.matches(np.load(fn))): make_synthetic.write(work, [seed], fam)
    z = np.load(fn); E, lab = z["E"][::stride].astype(float), z["lab"][::stride]
    return E, lab, xnl_ref.wrt(E)


def eps(P, Y, freq, ratio):
    df = freq * (ratio - 1 / ratio) / 2; p, y = P.sum(-1), Y.sum(-1)        # the direction step cancels
    return (np.abs(p - y) * df).sum(-1) / (np.abs(y) * df).sum(-1)


def table_check(work, E, Y):
    T = swrt.load_table(os.path.join(work, "table.npz")); full = swrt.Rows(T)
    P = np.concatenate([full.snl(E[s:s + 25]) for s in range(0, len(E), 25)])
    mr = np.abs(P - Y).max((1, 2)) / np.abs(Y).max((1, 2))
    r = dict(n=len(E), median=float(np.median(mr)), max=float(mr.max()))
    print(f"full table vs WRT on {len(E)} held-out spectra: max |P - Y| / max |Y| median {r['median']:.1e}, max {r['max']:.1e}", flush=True)
    return r


def operator_check(work, ratio, ops, E, lab, Y):
    T = swrt.load_table(os.path.join(work, "table.npz")); out = {}
    for op in ops:
        S = swrt.SWRT(op); P = np.concatenate([S.snl(E[s:s + 25]) for s in range(0, len(E), 25)]); e = eps(P, Y, T["freq"], ratio)
        out[os.path.basename(op)] = r = dict(K=S.K, eps_median=float(np.median(e)), eps_p95=float(np.percentile(e, 95)),
                                             eps_median_single_peak=float(np.median(e[lab == 0])), eps_median_two_peak=float(np.median(e[lab == 1])))
        print(f"{os.path.basename(op)}: eps median {r['eps_median']:.3f} (single peak {r['eps_median_single_peak']:.3f}, "
              f"two peaks {r['eps_median_two_peak']:.3f}), 95th percentile {r['eps_p95']:.3f}", flush=True)
    return out


def run(work, ratio, ops=(), seed=1000, stride=5, fam=None):
    E, lab, Y = heldout(work, seed, stride, fam)
    out = {"spectra": dict(seed=seed, stride=stride, n=len(E)), "table_vs_wrt_maxrel": table_check(work, E, Y),
           "operators": operator_check(work, ratio, ops, E, lab, Y)}
    json.dump(out, open(os.path.join(work, "check.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--f1", type=float, required=True); ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--nf", type=int, required=True); ap.add_argument("--ndir", type=int, required=True)
    ap.add_argument("--work", required=True); ap.add_argument("--ops", nargs="*", default=[])
    ap.add_argument("--seed", type=int, default=1000); ap.add_argument("--stride", type=int, default=5)
    a = ap.parse_args()
    freq, _ = xnl_ref.grid(a.f1, a.ratio, a.nf, a.ndir); xnl_ref.init(freq, a.ndir, os.path.join(a.work, "xnl"))
    run(a.work, a.ratio, a.ops, a.seed, a.stride)
