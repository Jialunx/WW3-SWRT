"""Make an SWRT operator for a WAVEWATCH III spectral grid in one go, by the steps of the paper (see README.md):
  1  WRT (Xnl) on the grid: its quadrature database and run-time constants                          xnl_ref.py
  2  the table of WRT's working rows, and its check against WRT on held-out practice spectra         build_table.py, check_operator.py
  3  practice spectra, 1500 per seed                                                                make_synthetic.py
  4  Gram matrices, greedy search with swaps, least-squares gains for each K                         build_swrt.py
  5  the file that NL6 reads                                                                         export_operator.py
  6  the operators' errors against WRT on the held-out spectra, and NL6 itself on them               check_operator.py, check_nl6.py
Grid: frequencies f1 * ratio^i (i = 0 .. nf-1), ndir directions over the full circle, deep water. Give the f1 (FREQ%F1) and
ratio (FREQ%XFR) of ww3_grid.nml; NL6 stops unless its grid matches the operator's. The work directory keeps every
intermediate file and is tied to one grid and one set of practice-spectrum ranges (grid.json): the practice spectra and the
Gram matrices are reused only when they belong to them.
Constants (--constants): 'paper' takes k, cg and dk from Xnl's functions at 1000 m, as the paper did, and reproduces its
files bit for bit; 'xnl' takes Xnl's own run-time arrays (at Xnl's deep-water depth of 2000 m), which agree better with WRT
when the grid starts at low frequencies; 'auto' (default) is 'paper' from 0.033 Hz up and 'xnl' below.
usage: python3 tools/make_operator.py --f1 0.03453 --ratio 1.1 --nf 40 --ndir 24 --work work_40x24 [--ks 30] [--seeds ...]
       ->  WORK/swrt_K30.ww3op: copy it to the run directory as swrt_operator.ww3op, or set WW3_SWRT_OPERATOR to its path"""
import os, sys, json, time, shutil, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import xnl_ref, build_table, make_synthetic, build_swrt, export_operator, check_operator, check_nl6

PAPER_GRID = (0.03453, 1.1, 40, 24)
PAPER_BQF = "01cfacad7cadec15eaeb0b0c21db9d43290f0d9d70fb26bdb4e95814d6a0d5c9"   # sha256 of the slots WRT reads (build_table.valid_hash)
F1_AUTO = 0.033


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--f1", type=float, required=True, help="lowest frequency (Hz), FREQ%%F1 of ww3_grid.nml")
    ap.add_argument("--ratio", type=float, required=True, help="frequency ratio, FREQ%%XFR")
    ap.add_argument("--nf", type=int, required=True, help="number of frequencies, FREQ%%NK")
    ap.add_argument("--ndir", type=int, required=True, help="number of directions, FREQ%%NTH")
    ap.add_argument("--work", required=True, help="work directory (one per grid)")
    ap.add_argument("--ks", default="30", help="numbers of kept configurations, e.g. 10,30 (default 30)")
    ap.add_argument("--seeds", default="11,1,2,3,4,5,6,7,8,9", help="practice-spectrum seeds (default: the paper's ten)")
    ap.add_argument("--fp-range", default="0.05,0.30", help="peak frequencies of the single-peak practice spectra (Hz)")
    ap.add_argument("--swell-range", default="0.05,0.12", help="peak frequencies of the swells of the two-peak spectra (Hz)")
    ap.add_argument("--constants", default="auto", choices=("auto", "paper", "xnl"), help="source of k, cg and dk (see above)")
    ap.add_argument("--no-nl6-check", action="store_true", help="skip the standalone check of W3SNL6 (needs gfortran)")
    a = ap.parse_args(); t0 = time.time(); work = os.path.abspath(a.work); os.makedirs(work, exist_ok=True)
    Ks = tuple(int(k) for k in a.ks.split(",")); seeds = [int(s) for s in a.seeds.split(",")]
    rng = lambda s: tuple(float(x) for x in s.split(","))
    if not a.no_nl6_check and not shutil.which("gfortran"): sys.exit("gfortran is needed for the NL6 check (or give --no-nl6-check)")
    const = a.constants if a.constants != "auto" else ("paper" if a.f1 >= F1_AUTO else "xnl")
    grid = dict(f1=a.f1, ratio=a.ratio, nf=a.nf, ndir=a.ndir, fp_range=list(rng(a.fp_range)), swell_range=list(rng(a.swell_range)), constants=const)
    gfile = os.path.join(work, "grid.json")
    if os.path.exists(gfile) and json.load(open(gfile)) != grid:
        sys.exit(f"{work} was made for {json.load(open(gfile))}; use another --work for {grid}")
    json.dump(grid, open(gfile, "w"), indent=1)

    print("1  WRT (Xnl) on the grid", flush=True)
    freq, _ = xnl_ref.grid(a.f1, a.ratio, a.nf, a.ndir); xnl_ref.init(freq, a.ndir, os.path.join(work, "xnl")); C = xnl_ref.consts(const)
    print(f"   {a.nf} frequencies {freq[0]:.5f}-{freq[-1]:.4f} Hz (ratio {a.ratio}), {a.ndir} directions; Xnl: {C['nlocus0']} points per locus, "
          f"no lumping, no Gauss-Legendre, filter (k ratio {C['qf_krat']}, {C['qf_dmax']} deg, density {C['qf_frac']:.2f}); constants '{const}'", flush=True)

    print("2  table of WRT's rows, and its check against WRT", flush=True)
    t, s = build_table.build(os.path.join(work, "xnl", "quad99999.bqf"), C); build_table.save(t, s, work)
    if (a.f1, a.ratio, a.nf, a.ndir) == PAPER_GRID and s["bqf_valid_sha256"] != PAPER_BQF:
        print("warning: Xnl's database differs from the paper's on its grid (compiler or flags, see build_xnl.sh); "
              "the operators will not match the paper's files bit for bit", flush=True)
    fam = make_synthetic.family(t, rng(a.fp_range), rng(a.swell_range))
    E, lab, Y = check_operator.heldout(work, fam=fam); tc = check_operator.table_check(work, E, Y)
    if tc["max"] > check_operator.TABLE_TOL: sys.exit(f"the table differs from WRT by {tc['max']:.1e} (more than {check_operator.TABLE_TOL:.0e}): stop")

    print("3  practice spectra", flush=True)
    todo = [sd for sd in seeds if not (os.path.exists(p := os.path.join(work, f"synthetic_s{sd}.npz")) and fam.matches(np.load(p)))]
    make_synthetic.write(work, todo, fam)

    print("4  selection and gains", flush=True)
    build_swrt.build(work, seeds, Ks, c_check=True)

    print("5  files for NL6", flush=True)
    ops = [os.path.join(work, f"swrt_K{K}.npz") for K in Ks]; ww3ops = [export_operator.export(op) for op in ops]

    print("6  errors against WRT on the held-out spectra, and NL6", flush=True)
    oc = check_operator.operator_check(work, a.ratio, ops, E, lab, Y)
    summary = {"grid": grid, "table": s, "table_vs_wrt_maxrel": tc, "operators": oc, "nl6_vs_python_maxrel": None, "seeds": seeds}
    out = os.path.join(work, "make_operator.json"); json.dump(dict(summary, minutes=(time.time() - t0) / 60), open(out, "w"), indent=1)
    if not a.no_nl6_check:
        summary["nl6_vs_python_maxrel"] = {os.path.basename(w): check_nl6.run(op, w, a.f1, a.ratio, a.nf, a.ndir, E[:50]) for op, w in zip(ops, ww3ops)}
        json.dump(dict(summary, minutes=(time.time() - t0) / 60), open(out, "w"), indent=1)
    print(f"done in {(time.time() - t0) / 60:.0f} min: " + ", ".join(ww3ops))


if __name__ == "__main__":
    main()
