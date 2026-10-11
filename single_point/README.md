# Single-point simulations

These are the six 48-h single-point simulations of the paper (its section 4 and Table 1). Each has one sea point at a depth of 1000 m in a 3 × 3 grid, a uniform wind and no propagation. The source terms are ST4 and LN1, and the spectral grid is that of the operators: 40 frequencies from 0.03453 Hz with a ratio of 1.1, and 24 directions. The spectrum and the source terms are written every 30 min.

| Directory | Scenario |
|---|---|
| `g20` | growth from calm under a wind of 20 m/s |
| `turn90` | the same wind, turning by 90° at 24 h |
| `decay20` | the same wind, dropping to 1 m/s at 24 h |
| `swell_cross` | swell (0.07 Hz, Hs 2 m, cos^20) with a wind of 10 m/s at 90° to it |
| `swell_opp` | swell (0.08 Hz, Hs 1.5 m, cos^20) with an opposing wind of 15 m/s |
| `swell_wide` | swell (0.06 Hz, Hs 3 m, cos^10) under a wind of 1 m/s |

Each directory holds the inputs of `ww3_grid`, `ww3_strt`, `ww3_shel` and `ww3_ounp`. They are the same for every S_nl option.

## Build

The single-point simulations use two settings that differ from the source in `ww3_src/`, which is set up for the global test:

- the semi-implicit integration of the source terms, with the diagonal term of S_nl, where `ww3_src/` integrates explicitly (`OFFSET = 0` in `w3srcemd.F90`);
- the DIA of WAVEWATCH III, where the DIA of `ww3_src/` skips spectra whose largest action density is below 1e-4.

`single_point.patch` restores both. Apply it in a separate copy of the archive, because the global test needs the unpatched source. In the top directory of that copy:

```
patch -p1 < single_point/single_point.patch
for n in 6 2 1; do                    # NL6 = SWRT, NL2 = WRT (Xnl), NL1 = DIA
  sed "s/NL6/NL$n/" single_point/switch > switch_nl$n
  cmake -S ww3_src -B build_nl$n -DSWITCH=$PWD/switch_nl$n -DCMAKE_BUILD_TYPE=Release
  make -C build_nl$n -j 8 ww3_grid ww3_strt ww3_shel ww3_ounp
done
```

`single_point/switch` is the switch file of the global test with `SHRD` (one process) in place of `DIST MPI`.

## Run

```
bash single_point/run.sh g20 build_nl2/bin                                    # WRT
bash single_point/run.sh g20 build_nl1/bin                                    # DIA
bash single_point/run.sh g20 build_nl6/bin                                    # SWRT, K = 30
bash single_point/run.sh g20 build_nl6/bin inputs/operators/swrt_K10.ww3op    # SWRT, K = 10
```

Each run has its own directory, `single_point/run_<scenario>_<build directory>`, with the operator name added for an operator given explicitly. The file `ww3.202501_src.nc` there holds the spectrum and the source terms every 30 min. The test spectra of the paper are spectra of the runs with WRT.

## Check

Built this way with gfortran 13.3 and run on one core, the runs with WRT and with DIA reproduce the point output of the paper's runs bit for bit in all six scenarios. With SWRT (K = 10, 30 and 100), the spectra differ from those of the paper's runs by less than 1e-6 of their largest value. The reason is that `ww3_src/` multiplies the transfer of SWRT by the WAM depth factor. At a depth of 1000 m this factor differs from 1 by less than 1e-5, and the paper's runs did not apply it.
