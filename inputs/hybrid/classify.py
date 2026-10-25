"""The rule of the hybrid in the paper's Discussion: the class of a spectrum E(f, theta) on the paper's grid (40 frequencies from 0.03453 Hz
with a ratio of 1.1, 24 equally spaced directions) and the configuration set of that class.
    sigma_w   circular directional spread of D(theta) = sum_f E(f, theta) f^4 df (the weight stresses the high-frequency waves)
    r         fraction of the energy more than 120 deg from the mean direction of D
    opposing seas (X): r >= 0.05;  narrow swell (S): otherwise sigma_w < 20 deg;  broad seas (W): all other spectra
Both quantities do not depend on the units of E or on the direction of the first bin.
usage: python3 classify.py spectra.npy   (E of shape (n, 40, 24) or (40, 24))  ->  one line per spectrum: sigma_w, r, class, operator file"""
import sys
import numpy as np

F = 0.03453 * 1.1 ** np.arange(40); DF = np.gradient(F); TH = np.deg2rad(np.arange(24) * 15.0)
SPREAD_DEG, ANGLE_DEG, FRAC = 20.0, 120.0, 0.05
SETS = {"S": "set_S.ww3op", "W": "set_W.ww3op", "X": "set_X.ww3op"}


def spread_w(E):
    D = (np.asarray(E, float) * (F ** 4 * DF)[None, :, None]).sum(1); D = D / D.sum(1, keepdims=True)
    R = np.hypot((D * np.cos(TH)).sum(1), (D * np.sin(TH)).sum(1)); return np.degrees(np.sqrt(2 * (1 - R)))


def opposed_frac(E, beyond=ANGLE_DEG):
    E = np.asarray(E, float); D = (E * (F ** 4 * DF)[None, :, None]).sum(1)
    tw = np.arctan2((D * np.sin(TH)).sum(1), (D * np.cos(TH)).sum(1))
    far = np.degrees(np.abs(np.angle(np.exp(1j * (TH[None, :] - tw[:, None]))))) > beyond
    e = (E * DF[None, :, None]).sum(1); return (e * far).sum(1) / e.sum(1)


def classes(E):
    E = np.asarray(E, float).reshape(-1, len(F), len(TH))
    return np.where(opposed_frac(E) >= FRAC, "X", np.where(spread_w(E) < SPREAD_DEG, "S", "W"))


if __name__ == "__main__":
    E = np.load(sys.argv[1]).astype(float).reshape(-1, len(F), len(TH))
    for s, r, c in zip(spread_w(E), opposed_frac(E), classes(E)): print(f"{s:7.2f} {r:7.4f} {c} {SETS[c]}")
