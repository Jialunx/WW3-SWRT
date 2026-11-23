"""Screen of candidate sets (candidates.py) in the six single-point training simulations of single_point/ (t_g12, t_g16, t_g25,
t_turn60, t_swell45, t_swell135), the step that retained SWRT9 among the candidates in the paper (its section 3d).

Each candidate runs 48 h in each training simulation with a single-point NL6 build; WRT (an NL2 build) is the reference, and
DIA (an NL1 build, optional) the basis of the step ratio. Per candidate and simulation, the screen records
  steps   the number of time steps of the source terms: the calls of the S_nl routine in ww3_shel, counted with gdb;
  error   the largest relative difference of H_s from that with WRT from 1 h on (spectra every 30 min).
The candidate with the fewest time steps in the six simulations together is retained (the paper's rule). The table also gives
the average over the six simulations of the largest error in H_s. As in the paper's screen, SWRT runs without the weak-wind
rule (WW3_SWRT_TOPM=0; --topm 2 applies it), and the paper's screen used the flush-to-zero library (--ftz build/ftz.so).
The builds are those of single_point/README.md (gfortran: gdb finds the routines by their gfortran names).
usage: python3 tools/screen.py --candidates work_40x24/candidates_K9/candidates.json --nl6 build_nl6/bin --nl2 build_nl2/bin
         [--nl1 build_nl1/bin] [--ftz build/ftz.so] [--jobs 4] [--keep-runs]
       ->  screen.json (all numbers) and screen.txt (the table) next to candidates.json"""
import os, re, sys, json, time, shutil, argparse, subprocess
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__)); SP = os.path.join(os.path.dirname(HERE), "single_point")
SIMS = ("t_g12", "t_g16", "t_g25", "t_turn60", "t_swell45", "t_swell135")
SNL = {"nl6": "__w3snl6md_MOD_w3snl6", "nl1": "__w3snl1md_MOD_w3snl1"}     # S_nl routines of SWRT and DIA (gfortran names)


def hs_series(nc):
    """time (h from the start) and H_s of the point output: 4 sqrt(sum E dtheta df), df = f (r - 1/r) / 2"""
    with Dataset(nc) as d:
        E = np.asarray(d["efth"][:, 0], float); f = np.asarray(d["frequency"][:], float); t = np.asarray(d["time"][:], float)
        nth = len(d["direction"][:])
    r = f[1] / f[0]; df = f * (r - 1 / r) / 2
    return (t - t[0]) * 24, 4 * np.sqrt((E.sum(-1) * (2 * np.pi / nth) * df).sum(-1))


def run(sim, bindir, kind, op, d, ftz, topm, keep):
    """one 48-h training simulation in a fresh copy d; returns the steps (None for WRT) and H_s every 30 min"""
    bindir, d = os.path.abspath(bindir), os.path.abspath(d)
    shutil.rmtree(d, ignore_errors=True); shutil.copytree(os.path.join(SP, sim), d)
    e = {k: v for k, v in os.environ.items() if not k.startswith("WW3_") and k != "LD_PRELOAD"}; e["OMP_NUM_THREADS"] = "1"
    if op: e.update(WW3_SWRT_OPERATOR=os.path.abspath(op), WW3_SWRT_TOPM=str(topm))
    if ftz: e["LD_PRELOAD"] = os.path.abspath(ftz)
    sh = lambda exe: subprocess.run([os.path.join(bindir, exe)], cwd=d, env=e, stdout=open(os.path.join(d, exe + ".log"), "w"),
                                    stderr=subprocess.STDOUT, check=True)
    sh("ww3_grid"); sh("ww3_strt"); steps = None
    if kind in SNL:
        p = subprocess.run(["gdb", "-batch", "-ex", "set pagination off", "-ex", f"break {SNL[kind]}", "-ex", "ignore 1 1000000000",
                            "-ex", "run > ww3_shel.log 2>&1", "-ex", "info breakpoints", os.path.join(bindir, "ww3_shel")],
                           cwd=d, env=e, capture_output=True, text=True)
        m = re.search(r"breakpoint already hit (\d+) time", p.stdout)
        if not m or "exited normally" not in p.stdout:
            raise RuntimeError(f"{d}: ww3_shel failed, or gdb did not find {SNL[kind]} (gfortran build?)\n{p.stdout[-800:]}")
        steps = int(m.group(1))
    else:
        sh("ww3_shel")
    sh("ww3_ounp"); t, hs = hs_series(os.path.join(d, "ww3.202501_src.nc"))
    if not keep: shutil.rmtree(d)
    return steps, t.tolist(), hs.tolist()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--candidates", required=True, help="candidates.json of candidates.py")
    ap.add_argument("--nl6", required=True, help="bin/ of the single-point NL6 build (SWRT)")
    ap.add_argument("--nl2", required=True, help="bin/ of the single-point NL2 build (WRT, the reference)")
    ap.add_argument("--nl1", default=None, help="bin/ of the single-point NL1 build (DIA, for the step ratio; optional)")
    ap.add_argument("--ftz", default=None, help="flush-to-zero library for LD_PRELOAD (build/ftz.so; the paper's screen used it)")
    ap.add_argument("--topm", type=int, default=0, help="WW3_SWRT_TOPM of the SWRT runs (default 0, as in the paper's screen)")
    ap.add_argument("--jobs", type=int, default=4, help="simulations at a time (default 4)")
    ap.add_argument("--keep-runs", action="store_true", help="keep the run directories")
    a = ap.parse_args(); t0 = time.time()
    if shutil.which("gdb") is None: sys.exit("gdb is needed to count the time steps")
    cj = os.path.abspath(a.candidates); out = os.path.dirname(cj); C = json.load(open(cj))["candidates"]; rd = os.path.join(out, "screen_runs")
    jobs = [("WRT", s, a.nl2, "nl2", None) for s in SIMS] + ([("DIA", s, a.nl1, "nl1", None) for s in SIMS] if a.nl1 else [])
    jobs += [(c["name"], s, a.nl6, "nl6", os.path.join(out, c["ww3op"])) for c in C for s in SIMS]
    R = {}
    with ProcessPoolExecutor(a.jobs) as ex:
        fut = {ex.submit(run, s, os.path.abspath(b), k, op, os.path.join(rd, f"{n}_{s}"), a.ftz, a.topm, a.keep_runs): (n, s)
               for n, s, b, k, op in jobs}
        for f, (n, s) in fut.items(): R.setdefault(n, {})[s] = f.result()
    ref = {s: np.array(R["WRT"][s][2]) for s in SIMS}; m1 = {s: np.array(R["WRT"][s][1]) >= 1 - 1e-6 for s in SIMS}
    dia = sum(R["DIA"][s][0] for s in SIMS) if a.nl1 else None; res = {}
    for c in C:
        r = R[c["name"]]; err = {s: float(np.max(np.abs(np.array(r[s][2])[m1[s]] / ref[s][m1[s]] - 1))) for s in SIMS}
        st = {s: r[s][0] for s in SIMS}; tot = sum(st.values())
        res[c["name"]] = dict(e2=c["e2"], offsets=c["offsets"], start=c["start"], limit=c["limit"], seed=c["seed"], steps=st,
                              steps_total=tot, steps_rel_dia=tot / dia if dia else None, hs_error=err,
                              hs_error_mean=float(np.mean(list(err.values()))), hs_error_max=float(max(err.values())))
    keep = min(res, key=lambda n: (res[n]["steps_total"], res[n]["hs_error_mean"]))
    json.dump(dict(simulations=SIMS, rule="fewest time steps in the six training simulations together", retained=keep,
                   dia_steps=({s: R["DIA"][s][0] for s in SIMS} if a.nl1 else None), topm=a.topm, ftz=bool(a.ftz),
                   candidates=res, wrt_hs={s: R["WRT"][s][2] for s in SIMS}, minutes=(time.time() - t0) / 60),
              open(os.path.join(out, "screen.json"), "w"), indent=1)
    lines = [f"{'':5s} {'e2':>7s} {'L':>3s} {'steps':>6s} " + ("x DIA " if dia else "") + f"{'mean':>6s} {'max':>6s}  largest H_s error per simulation (%)"]
    for n in sorted(res, key=lambda n: (res[n]["steps_total"], res[n]["hs_error_mean"])):
        v = res[n]
        lines.append(f"{n:5s} {v['e2']:7.4f} {v['offsets']:3d} {v['steps_total']:6d} " + (f"{v['steps_rel_dia']:5.2f} " if dia else "")
                     + f"{100 * v['hs_error_mean']:5.2f}% {100 * v['hs_error_max']:5.2f}%  "
                     + " ".join(f"{100 * v['hs_error'][s]:5.2f}" for s in SIMS) + ("   <- retained" if n == keep else ""))
    open(os.path.join(out, "screen.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
    print(f"retained: {keep} -> {os.path.join(out, C[[c['name'] for c in C].index(keep)]['ww3op'])} ({(time.time() - t0) / 60:.0f} min)")


if __name__ == "__main__":
    main()
