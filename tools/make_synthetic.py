"""Synthetic practice spectra for the selection (build_swrt.py): a JONSWAP family generated from formulas and one random seed
on the grid of the table; no model data and no WRT answers. Per seed (the paper pools seeds 11 and 1-9, 15 000 spectra):
  900 single-peak spectra: fp log-uniform 0.05-0.30 Hz, gamma uniform 1-7, directional spread at the peak log-uniform 4-60 deg
      (cos^2s(dtheta/2) with s solved on the direction grid), random mean direction, Hs 2 m;
  600 two-peak spectra: swell fp1 log-uniform 0.05-0.12 Hz, gamma 1-7, spread 4-20 deg; wind sea fp2 = (1.3-3.0) fp1, gamma 1-5,
      spread 20-50 deg; independent random directions; swell scaled by 0.2-5 relative to the wind sea (each Hs 1 m before scaling).
The frequency ranges are those of the paper (grid 0.0345-1.42 Hz). For another frequency range, set them with --fp-range and
--swell-range so that the peaks lie well inside the grid (the paper's lowest peak, 0.05 Hz, is four bins above the lowest
frequency, and the highest wind-sea peak, 0.36 Hz, is 14 bins below the highest); a warning is printed otherwise.
Each file also holds the grid and the ranges, so that a file made for another grid or other ranges is not reused.
usage: python3 make_synthetic.py --work DIR [--seeds 11,1,2,3,4,5,6,7,8,9]  ->  DIR/synthetic_s{seed}.npz (needs DIR/table.npz)"""
import os, argparse
import numpy as np

G0 = 9.81


class Family:
    def __init__(self, freq, nd, fp_range=(0.05, 0.30), swell_range=(0.05, 0.12)):
        self.f = np.asarray(freq, float); self.df = np.gradient(self.f); self.nd = nd
        self.th = np.deg2rad(np.arange(nd) * (360.0 / nd)); self.fp_range, self.swell_range = tuple(fp_range), tuple(swell_range)
        lo, hi = min(fp_range[0], swell_range[0]), max(fp_range[1], 3.0 * swell_range[1])
        if lo < self.f[min(3, len(self.f) - 1)] or hi > self.f[max(len(self.f) - 9, 0)]:
            print(f"warning: practice-spectrum peaks from {lo:.3f} to {hi:.3f} Hz are close to or outside the grid "
                  f"({self.f[0]:.4f}-{self.f[-1]:.3f} Hz); set --fp-range and --swell-range", flush=True)

    def matches(self, z):
        """whether a stored practice-spectrum file belongs to this grid and these ranges"""
        return ("freq" in z and z["E"].shape[1:] == (len(self.f), self.nd) and np.array_equal(z["freq"], self.f)
                and tuple(z["fp_range"]) == self.fp_range and tuple(z["swell_range"]) == self.swell_range)

    def spread_deg(self, D):
        D = D / D.sum(); C, S = (D * np.cos(self.th)).sum(), (D * np.sin(self.th)).sum(); return np.degrees(np.sqrt(2 * (1 - np.hypot(C, S))))

    def dirdist(self, theta_p, target):
        """cos^(2s)((theta - theta_p)/2) on the direction grid, s solved so that the discrete spread equals target (deg)"""
        base = np.abs(np.cos(0.5 * (self.th - np.deg2rad(theta_p)))); lo, hi = 0.5, 4000.0     # |cos| = cos of the wrapped half-difference
        for _ in range(60):
            s = np.sqrt(lo * hi); D = base ** (2 * s)
            lo, hi = (s, hi) if self.spread_deg(D + 1e-300) > target else (lo, s)
        D = base ** (2 * s); return D / (D.sum() * 2 * np.pi / self.nd)

    def jonswap1d(self, fp, gamma, Hs):
        f = self.f; sig = np.where(f <= fp, 0.07, 0.09); rr = np.exp(-((f - fp) ** 2) / (2 * (sig * fp) ** 2))
        S1 = G0 ** 2 / ((2 * np.pi) ** 4 * f ** 5) * np.exp(-1.25 * (fp / f) ** 4) * gamma ** rr
        return S1 * (Hs / 4) ** 2 / (S1 * self.df).sum()

    def make(self, seed):
        rng = np.random.default_rng(seed); E, lab = [], []; (a, b), (c, d) = np.log(self.fp_range), np.log(self.swell_range)
        for _ in range(900):                                   # single peak
            fp = np.exp(rng.uniform(a, b)); gam = rng.uniform(1.0, 7.0); sp = np.exp(rng.uniform(np.log(4.0), np.log(60.0)))
            E.append(self.jonswap1d(fp, gam, 2.0)[:, None] * self.dirdist(rng.uniform(0, 360), sp)[None, :]); lab.append(0)
        for _ in range(600):                                   # swell + wind sea
            fp1 = np.exp(rng.uniform(c, d)); fp2 = fp1 * rng.uniform(1.3, 3.0)
            e1 = self.jonswap1d(fp1, rng.uniform(1, 7), 1.0)[:, None] * self.dirdist(rng.uniform(0, 360), np.exp(rng.uniform(np.log(4.0), np.log(20.0))))[None, :]
            e2 = self.jonswap1d(fp2, rng.uniform(1, 5), 1.0)[:, None] * self.dirdist(rng.uniform(0, 360), rng.uniform(20.0, 50.0))[None, :]
            E.append(e1 * rng.uniform(0.2, 5.0) + e2); lab.append(1)
        return np.stack(E).astype(np.float32), np.array(lab)


def family(T, fp_range=(0.05, 0.30), swell_range=(0.05, 0.12)):
    from swrt import n_dirs
    return Family(T["freq"], n_dirs(T), fp_range, swell_range)


def write(work, seeds, fam):
    for seed in seeds:
        out = os.path.join(work, f"synthetic_s{seed}.npz")
        E, lab = fam.make(seed)
        np.savez(out, E=E, lab=lab, seed=seed, freq=fam.f, fp_range=fam.fp_range, swell_range=fam.swell_range, note=__doc__.split("\n")[0])
        print(f"seed {seed}: {len(E)} spectra ({(lab == 0).sum()} single-peak, {(lab == 1).sum()} two-peak) -> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", required=True); ap.add_argument("--seeds", default="11,1,2,3,4,5,6,7,8,9")
    ap.add_argument("--fp-range", default="0.05,0.30"); ap.add_argument("--swell-range", default="0.05,0.12")
    a = ap.parse_args(); T = dict(np.load(os.path.join(a.work, "table.npz")))
    rng_ = lambda s: tuple(float(x) for x in s.split(","))
    write(a.work, [int(s) for s in a.seeds.split(",")], family(T, rng_(a.fp_range), rng_(a.swell_range)))
