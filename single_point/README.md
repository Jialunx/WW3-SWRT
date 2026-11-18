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

The single-point simulations use the DIA of WAVEWATCH III, where the DIA of `ww3_src/`, set up for the global simulation, skips spectra whose largest action density is below 1e-4. `single_point.patch` restores it. (Up to v1.4, the patch also switched `ww3_src/` from explicit to semi-implicit integration; since v1.5 `ww3_src/` integrates semi-implicitly.) Apply it in a separate copy of the archive, because the global simulation needs the unpatched source. In the top directory of that copy:

```
patch -p1 < single_point/single_point.patch
for n in 6 2 1; do                    # NL6 = SWRT, NL2 = WRT (Xnl), NL1 = DIA
  sed "s/NL6/NL$n/" single_point/switch > switch_nl$n
  cmake -S ww3_src -B build_nl$n -DSWITCH=$PWD/switch_nl$n -DCMAKE_BUILD_TYPE=Release
  make -C build_nl$n -j 8 ww3_grid ww3_strt ww3_shel ww3_ounp
done
```

`single_point/switch` is the switch file of the global simulation with `SHRD` (one process) in place of `DIST MPI`.

## Run

```
bash single_point/run.sh g20 build_nl2/bin                                    # WRT
bash single_point/run.sh g20 build_nl1/bin                                    # DIA
bash single_point/run.sh g20 build_nl6/bin                                    # SWRT36 (inputs/swrt_K36.ww3op)
bash single_point/run.sh g20 build_nl6/bin inputs/swrt9.ww3op                 # SWRT9
bash single_point/run.sh g20 build_nl6/bin path/to/other.ww3op                # SWRT with another operator
```

`run.sh` sets the weak-wind rule of the paper (`WW3_SWRT_TOPM=2 WW3_SWRT_TOPM_UST=1`) unless these variables are already set; it acts only in NL6 builds. Each run has its own directory, `single_point/run_<scenario>_<build directory>`, with the operator name added for an operator given explicitly. The file `ww3.202501_src.nc` there holds the spectrum and the source terms every 30 min. The test spectra of the paper are spectra of the runs with WRT.

## Check

The check below was made with v1.3, before the weak-wind rule. Built this way with gfortran 13.3 and run on one core, the runs with WRT and with DIA reproduce the point output of the paper's runs bit for bit in all six scenarios (checked with v1.2; v1.3 changes only the code of SWRT). With SWRT (v1.3, K = 36), four of the six scenarios reproduce the paper's runs bit for bit, and in `swell_cross` and `swell_wide` the spectra differ from those of the paper's runs by less than 1e-6 of their largest value. The reason is that `ww3_src/` multiplies the transfer of SWRT by the WAM depth factor. At a depth of 1000 m this factor differs from 1 by up to 2.4e-5 in `swell_wide`, 2e-7 in `swell_cross` and less than 2e-9 in the other four scenarios, where it rounds to 1 in single precision. The paper's single-point runs did not apply it.
