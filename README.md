# WW3-SWRT

WAVEWATCH III (version 7.14) with the sparse Webb–Resio–Tracy (SWRT) source term for the nonlinear four-wave interaction, selected by the switch `NL6`.

SWRT retains K mirror-image pairs of configurations of the WRT integration space, each with an optimized gain. The operator of the paper, SWRT36, has K = 36 configurations, whose rows use 18 offsets of k3 from k1 (`inputs/swrt_K36.ww3op`). The faster SWRT9 of the paper has K = 9 configurations on 11 offsets (`inputs/swrt9.ww3op`). This version includes

- the kernel that groups the rows by resonance curve, used by default: the rows of one offset share N1 and N3, so that two sums per curve give the transfer and the diagonal term (set `WW3_SWRT_KERNEL=vector` for the vectorized row loop of v1.2, or `WW3_SWRT_KERNEL=reference` for the original loop),
- S_nl at all frequencies of the grid, as with DIA and WRT (`NKLIM = 0` in `w3srcemd.F90`),
- the WAM depth scaling in finite depth, as Xnl applies it to WRT with `IQTYPE = 2` (the same function of kd as for DIA, `NL1`),
- the diagonal term for implicit time stepping (set `WW3_SWRT_DIAG=0` to switch it off),
- the semi-implicit integration of the source terms (`OFFSET = 1` in `w3srcemd.F90`) with a dynamic time step of 15–900 s, as in all simulations of the paper,
- a weak-wind rule for the two highest prognostic frequencies (on by default, as in all SWRT runs of the paper; `WW3_SWRT_TOPM=0` switches it off): where ST4 sets the highest prognostic frequency from the friction velocity, the diagonal term of SWRT at the two highest prognostic frequencies is raised to the value at which the semi-implicit scheme keeps the change of the spectrum in one step below the allowed change, so that these two frequencies no longer limit the time step.


## Contents

| Path | Content |
|---|---|
| `ww3_src/` | WAVEWATCH III source; SWRT is in `model/src/w3snl6md.F90` |
| `switch` | WAVEWATCH III switches (`NL6` selects SWRT) |
| `inputs/swrt_K36.ww3op` | SWRT operator of the paper, SWRT36 (K = 36) |
| `inputs/swrt9.ww3op` | the faster operator of the paper, SWRT9 (K = 9 on 11 offsets) |
| `inputs/ftz.c` | flush-to-zero library for `LD_PRELOAD` |
| `inputs/Global-real.depth`, `inputs/points.list` | bathymetry and output points of the 1° global simulation of the paper |
| `namelists/` | namelists of the global simulation |
| `scripts/` | build, run-directory and Slurm scripts, written for an HPC cluster with environment modules (adapt module names and paths) |
| `tools/` | scripts that make the SWRT operator for another spectral grid, with the search of the paper (`tools/README.md`) |
| `single_point/` | setups of the six single-point simulations of the paper, with a patch for their build (`single_point/README.md`) |

The operators of v1.2 (K = 30, and those for other K in `inputs/operators/`) came from an earlier selection and are not part of this version. They remain in the v1.2 archive. The configuration sets of the hybrid (v1.3, v1.4) are not part of this version; they remain in the v1.4 archive.

Two inputs of the global simulation are not included because of their size: the initial condition `restart_ic_20250101.ww3` (138 MB) and the ERA5 wind forcing `wind_1deg_global_20250101to15.nc` (351 MB, 1°, 1–15 January 2025). Prepare the wind file from ERA5 10-m winds, and provide an initial condition or start from rest, before running `scripts/02_prep_run.sh`.

## Spectral grid

The operator file belongs to one spectral grid: 40 frequencies from 0.03453 Hz with a ratio of 1.1, and 24 directions. The model checks this grid when it reads the operator. For another grid, make the operator with `tools/make_operator.py`, which follows the steps of the paper: WRT's quadrature for the grid, the table of its rows, practice spectra, and the selection of the configurations, optionally under a limit on their offsets, and of their gains (see `tools/README.md`).

## Build

```
bash scripts/01_build.sh        # CMake with -DSWITCH=switch; also builds build/ftz.so
```

## Run

```
bash scripts/02_prep_run.sh                 # creates the run directory with SWRT36
bash scripts/02_prep_run.sh inputs/swrt9.ww3op   # ... or with SWRT9
export WW3_SWRT_OPERATOR=/path/to/swrt_operator.ww3op   # default: swrt_operator.ww3op in the run directory
LD_PRELOAD=/path/to/build/ftz.so mpirun -n <N> -x WW3_SWRT_OPERATOR /path/to/build/bin/ww3_shel   # LD_PRELOAD is optional
```

## License

WAVEWATCH III is distributed under the GNU Lesser General Public License (see `LICENSE.md`). The SWRT changes are distributed under the same terms.

## Reference

Chen, J., T. A. A. Adcock, Q. Liu, and T. Tang: A new nonlinear four-wave interaction scheme for ocean wave models: Webb–Resio–Tracy accuracy at operational speed. Manuscript in preparation.
