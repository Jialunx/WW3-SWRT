# WW3-SWRT

WAVEWATCH III (version 7.14) with the sparse Webb–Resio–Tracy (SWRT) source term for the nonlinear four-wave interaction, selected by the switch `NL6`.

SWRT keeps K mirror-image pairs of configurations of the WRT integration space, each with an optimized gain. The operator here has K = 30. This version includes

- the vectorized loop over directions, used by default (set `WW3_SWRT_KERNEL=reference` for the original loop),
- the WAM depth scaling in finite depth, the same factor as for DIA (`NL1`) and for WRT with `IQTYPE = 2`,
- the diagonal term for implicit time stepping (set `WW3_SWRT_DIAG=0` to switch it off),
- an optional library, `inputs/ftz.c`, that flushes subnormal floating-point numbers to zero when it is preloaded.

## Contents

| Path | Content |
|---|---|
| `ww3_src/` | WAVEWATCH III source; SWRT is in `model/src/w3snl6md.F90` |
| `switch` | WAVEWATCH III switches (`NL6` selects SWRT) |
| `inputs/swrt_K30.ww3op` | SWRT operator, K = 30 |
| `inputs/ftz.c` | flush-to-zero library for `LD_PRELOAD` |
| `inputs/Global-real.depth`, `inputs/points.list` | bathymetry and output points of a 1° global test |
| `namelists/` | namelists of the global test |
| `scripts/` | build, run-directory and Slurm scripts, written for an HPC cluster with environment modules (adapt module names and paths) |
| `tools/` | scripts that make the SWRT operator for another spectral grid (`tools/README.md`) |

Two inputs of the global test are not included because of their size: the initial condition `restart_ic_20250101.ww3` (138 MB) and the ERA5 wind forcing `wind_1deg_global_20250101to15.nc` (351 MB, 1°, 1–15 January 2025). Prepare the wind file from ERA5 10-m winds, and provide an initial condition or start from rest, before running `scripts/02_prep_run.sh`.

## Spectral grid

The operator file belongs to one spectral grid: 40 frequencies from 0.03453 Hz with a ratio of 1.1, and 24 directions. The model checks this grid when it reads the operator. For another grid, make the operator with `tools/make_operator.py`, which follows the steps of the paper: WRT's quadrature for the grid, the table of its rows, practice spectra, and the selection of the configurations and their gains (see `tools/README.md`).

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

Chen, J., T. A. A. Adcock, Q. Liu, R. Clark, and T. Tang: A sparse Webb–Resio–Tracy scheme of nonlinear energy transfer among wind waves for wave forecasting models. Manuscript in preparation.
