"""Standalone check of WAVEWATCH III's NL6 (ww3_src/model/src/w3snl6md.F90) with an operator file on its grid, before a model run.

The module is compiled with gfortran and WAVEWATCH III's GNU flags (release -O3), with the model's constants.F90 and minimal
stand-ins of the data modules (nl6_check/nl6_stubs.F90), and driven as W3SRCE drives it (nl6_check/nl6_driver.F90: SIG from
FR1 and XFR as ww3_grid builds it, A = E * CONX, S(E) = S / CONX, depth 10 km so that the WAM depth factor of W3SNL6 is 1).
On held-out practice spectra of the grid, S and the diagonal term D are compared with swrt.SWRT(operator.npz):
max |S_nl6 - S| / max |S| per spectrum must stay below 1e-5 (on the paper's grid it is below 1e-6).
usage: python3 check_nl6.py --op DIR/swrt_K30.npz [--ww3op DIR/swrt_K30.ww3op] --f1 0.03453 --ratio 1.1 --nf 40 --ndir 24
       (FR1 and XFR as in ww3_grid.nml: FREQ%F1 and FREQ%XFR; uses DIR/synthetic_s1000.npz, made by check_operator.py)"""
import os, sys, shutil, argparse, subprocess, tempfile
import numpy as np
import swrt

HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(os.path.dirname(HERE), "ww3_src", "model", "src")
FLAGS = "-cpp -DW3_NL6 -O3 -g -fno-second-underscore -ffree-line-length-none -ffpe-summary=none".split()   # WW3 GNU flags, release
TOL = 1e-5


def build(d):
    if not shutil.which("gfortran"): sys.exit("gfortran is needed for the NL6 check")
    srcs = [os.path.join(SRC, "constants.F90"), os.path.join(HERE, "nl6_check", "nl6_stubs.F90"), os.path.join(SRC, "w3snl6md.F90"),
            os.path.join(HERE, "nl6_check", "nl6_driver.F90")]
    objs = []
    for s in srcs:
        o = os.path.join(d, os.path.basename(s) + ".o"); objs.append(o)
        r = subprocess.run(["gfortran", *FLAGS, "-I", SRC, "-J", d, "-c", s, "-o", o], capture_output=True, text=True)
        if r.returncode: sys.exit(f"compile failed: {s}\n{r.stderr}")
    exe = os.path.join(d, "nl6_driver"); subprocess.run(["gfortran", *objs, "-o", exe], check=True)
    return exe


def run(op, ww3op, f1, ratio, nf, ndir, E):
    with tempfile.TemporaryDirectory() as d:
        exe = build(d); fin, fout = os.path.join(d, "in.bin"), os.path.join(d, "out.bin")
        np.ascontiguousarray(E, dtype="<f4").tofile(fin)
        env = dict(os.environ, WW3_SWRT_OPERATOR=os.path.abspath(ww3op)); env.pop("WW3_SWRT_DIAG", None)
        r = subprocess.run([exe, fin, fout, str(len(E)), repr(f1), repr(ratio), str(nf), str(ndir)], capture_output=True, text=True, env=env, cwd=d)
        if r.returncode: sys.exit(f"W3SNL6 stopped:\n{r.stdout}{r.stderr}")
        SD = np.fromfile(fout, "<f4").reshape(2, *E.shape).astype(float)
    S = swrt.SWRT(op); P, D = S.snl(E), S.diag(E)
    rel = lambda X, Y: np.abs(X - Y).max((1, 2)) / np.abs(Y).max((1, 2))
    rs, rd = rel(SD[0], P), rel(SD[1], D)
    print(f"W3SNL6 vs swrt.SWRT on {len(E)} spectra: S max rel {rs.max():.1e}, D max rel {rd.max():.1e} (tolerance {TOL:.0e})")
    if max(rs.max(), rd.max()) > TOL: sys.exit("NL6 check FAILED")
    return float(rs.max()), float(rd.max())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--op", required=True); ap.add_argument("--ww3op", default=None)
    ap.add_argument("--f1", type=float, required=True); ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--nf", type=int, required=True); ap.add_argument("--ndir", type=int, required=True)
    ap.add_argument("--n", type=int, default=50, help="number of held-out spectra")
    a = ap.parse_args()
    E = np.load(os.path.join(os.path.dirname(os.path.abspath(a.op)), "synthetic_s1000.npz"))["E"][:a.n].astype(float)
    run(a.op, a.ww3op or os.path.splitext(a.op)[0] + ".ww3op", a.f1, a.ratio, a.nf, a.ndir, E)
