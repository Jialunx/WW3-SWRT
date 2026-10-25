# WW3-SWRT

WAVEWATCH III (version 7.14) with the sparse Webb–Resio–Tracy (SWRT) source term for the nonlinear four-wave interaction, selected by the switch `NL6`.

SWRT retains K mirror-image pairs of configurations of the WRT integration space, each with an optimized gain. The operator of the paper has K = 36 configurations, whose rows use 18 offsets of k3 from k1 (`inputs/swrt_K36.ww3op`). This version includes

- the kernel that groups the rows by resonance curve, used by default: the rows of one offset share N1 and N3, so that two sums per curve give the transfer and the diagonal term (set `WW3_SWRT_KERNEL=vector` for the vectorized row loop of v1.2, or `WW3_SWRT_KERNEL=reference` for the original loop),
- S_nl at all frequencies of the grid, as with DIA and WRT (`NKLIM = 0` in `w3srcemd.F90`),
- the WAM depth scaling in finite depth, as Xnl applies it to WRT with `IQTYPE = 2` (the same function of kd as for DIA, `NL1`),
- the diagonal term for implicit time stepping (set `WW3_SWRT_DIAG=0` to switch it off),
- an optional library, `inputs/ftz.c`, that flushes subnormal floating-point numbers to zero when it is preloaded.

The code of `w3snl6md.F90` and `w3srcemd.F90` is that of the paper's global simulation and timings; v1.3 corrects only comments in `w3snl6md.F90` (the compiled code is identical).

## Contents

| Path | Content |
|---|---|
| `ww3_src/` | WAVEWATCH III source; SWRT is in `model/src/w3snl6md.F90` |
| `switch` | WAVEWATCH III switches (`NL6` selects SWRT) |
| `inputs/swrt_K36.ww3op` | SWRT operator of the paper, K = 36 |
| `inputs/hybrid/` | the three configuration sets of the hybrid in the paper's Discussion and the rule that assigns a spectrum to one of them (`inputs/hybrid/README.md`) |
| `inputs/ftz.c` | flush-to-zero library for `LD_PRELOAD` |
| `inputs/Global-real.depth`, `inputs/points.list` | bathymetry and output points of the 1° global simulation of the paper |
| `namelists/` | namelists of the global simulation |
| `scripts/` | build, run-directory and Slurm scripts, written for an HPC cluster with environment modules (adapt module names and paths) |
| `tools/` | scripts that make the SWRT operator for another spectral grid, with the search of the paper (`tools/README.md`) |
| `single_point/` | setups of the six single-point simulations of the paper, with a patch for their build (`single_point/README.md`) |

The operators of v1.2 (K = 30, and those for other K in `inputs/operators/`) came from an earlier selection and are not part of this version. They remain in the v1.2 archive.

Two inputs of the global simulation are not included because of their size: the initial condition `restart_ic_20250101.ww3` (138 MB) and the ERA5 wind forcing `wind_1deg_global_20250101to15.nc` (351 MB, 1°, 1–15 January 2025). Prepare the wind file from ERA5 10-m winds, and provide an initial condition or start from rest, before running `scripts/02_prep_run.sh`.

## Spectral grid

The operator file belongs to one spectral grid: 40 frequencies from 0.03453 Hz with a ratio of 1.1, and 24 directions. The model checks this grid when it reads the operator. For another grid, make the operator with `tools/make_operator.py`, which follows the steps of the paper: WRT's quadrature for the grid, the table of its rows, practice spectra, and the selection of the configurations, optionally under a limit on their offsets, and of their gains (see `tools/README.md`).

## Build

```
bash scripts/01_build.sh        # CMake with -DSWITCH=switch; also builds build/ftz.so
```

## Run

```
bash scripts/02_prep_run.sh     # creates the run directory
export WW3_SWRT_OPERATOR=/path/to/swrt_operator.ww3op   # default: swrt_operator.ww3op in the run directory
LD_PRELOAD=/path/to/build/ftz.so mpirun -n <N> /path/to/build/bin/ww3_shel   # LD_PRELOAD is optional
```

## License

WAVEWATCH III is distributed under the GNU Lesser General Public License (see `LICENSE.md`). The SWRT changes are distributed under the same terms.

## Reference

Chen, J., T. A. A. Adcock, Q. Liu, R. Clark, and T. Tang: A new nonlinear four-wave interaction scheme for ocean wave models: Webb–Resio–Tracy accuracy near DIA cost. Manuscript in preparation.
