# Making an SWRT operator for another spectral grid

The operator file `inputs/swrt_K30.ww3op` belongs to one spectral grid (40 frequencies from 0.03453 Hz with a ratio of 1.1,
24 directions), and NL6 stops if the model's grid differs. The scripts here make the operator for any grid of the same kind
by the steps of the paper:

1. WRT (Van Vledder's Xnl, the `NL2` of WAVEWATCH III) is set up on the grid. It writes its quadrature database
   `quad99999.bqf` and gives its run-time constants (`xnl_ref.py`).
2. The database is folded into the table of the rows that WRT uses at run time (`build_table.py`). Summed over all rows,
   the table reproduces WRT, and this is checked on held-out practice spectra (`check_operator.py`).
3. Practice spectra are generated from formulas, a JONSWAP family with one and two peaks, 1500 per seed (`make_synthetic.py`).
   No model output and no WRT result is used.
4. The K configurations (mirror-image pairs of rows) are chosen by a greedy search with swaps, and their gains by linear
   least squares, so that they reproduce the full table on the practice spectra (`build_swrt.py`).
5. The operator is written in the plain-text format that NL6 reads (`export_operator.py`).
6. Checks on held-out practice spectra: the error of each operator against WRT, and NL6 itself (`w3snl6md.F90`, compiled
   standalone) against the Python operator (`check_operator.py`, `check_nl6.py`).

## Requirements

Python 3 with NumPy, gfortran, and for Python 3.12 or later meson and ninja (`pip install meson ninja`), which `f2py` needs.
The checks below were made with gfortran 13.3, Python 3.12.3 and NumPy 1.26.4 on x86-64 Linux. The memory needed grows with
the grid: about 1.5 GB for the paper's grid and 3 GB for 36 frequencies and 36 directions.

## Use

```
bash tools/build_xnl.sh                       # once: Xnl as a Python module (tools/build/xnl_wrt*.so)
python3 tools/make_operator.py --f1 0.0345 --ratio 1.1 --nf 36 --ndir 36 --work work_36x36 --ks 30
```

`--f1`, `--ratio`, `--nf` and `--ndir` are `FREQ%F1`, `FREQ%XFR`, `FREQ%NK` and `FREQ%NTH` of `ww3_grid.nml`. The result is
`work_36x36/swrt_K30.ww3op`. Copy it to the run directory as `swrt_operator.ww3op`, or set `WW3_SWRT_OPERATOR` to its path.
Several values of K can be made at once (`--ks 10,30`).

The work directory keeps every intermediate file, and `make_operator.json` summarizes the table, the checks and the errors.
A work directory belongs to one grid and one set of practice-spectrum ranges (`grid.json`), and `make_operator.py` stops if
it is given another; use one work directory per grid. Within it, the practice spectra and the Gram matrices are reused
only if they were made for the same grid and table, so that other values of K take only the search, a few seconds.

The steps can also be run one by one. The scripts with a command line (`build_table.py`, `make_synthetic.py`,
`build_swrt.py`, `export_operator.py`, `check_operator.py`, `check_nl6.py`, `make_operator.py`) describe their arguments
with `-h`.

| Script | Step |
|---|---|
| `build_xnl.sh` | builds Xnl from `ww3_src/model/src` (`mod_constants.f90`, `mod_fileio.f90`, `serv_xnl4v5.f90`, `mod_xnl4v5.f90`) |
| `xnl_ref.py` | WRT on a grid: database, constants, transfer (module) |
| `bqf_decode.py` | reader of the database, in either byte order (module) |
| `build_table.py` | table of WRT's rows |
| `make_synthetic.py` | practice spectra |
| `swrt.py` | the rows and SWRT in Python, transfer and diagonal term (module) |
| `build_swrt.py` | selection and gains |
| `export_operator.py` | operator file for NL6 |
| `check_operator.py` | table and operators against WRT |
| `check_nl6.py`, `nl6_check/` | NL6 against the Python operator |
| `make_operator.py` | all steps |

## Run time

The Gram matrices of the 15 000 practice spectra take most of the time. For the paper's grid they took 17 minutes on an
otherwise idle Intel i9-11900K. The test runs below, made while several other jobs shared that machine, took 50 minutes for
the paper's grid and 68 minutes for 36 frequencies and 36 directions.

## Checks

- **The paper's operator is reproduced.** For the paper's grid, `make_operator.py` gives the same table (bit for bit),
  the same practice spectra, the same configurations and gains, and an operator file identical to `inputs/swrt_K30.ww3op`.
  The same holds for K = 10, 20 and 50 against the files of the paper's runs (not included here): the configurations are
  those of the paper, the gains agree to 1e-14, and the files are identical. The full table reproduces WRT to 5e-5 of the
  largest transfer on 300 held-out spectra, and NL6 reproduces the Python operator to 1e-6.
- **Another grid.** For 36 frequencies from 0.03453 Hz (ratio 1.1) and 36 directions, the table has 807 configurations
  (592 on the paper's grid), and the full table reproduces WRT to 4e-5. On held-out practice spectra, the median error of
  SWRT against WRT (`check_operator.py`) is 0.08 with K = 30 and 0.36 with K = 10, against 0.06 and 0.20 on the paper's
  grid. A grid with more directions thus needs more configurations for the same error. NL6 reproduces the Python operator
  to 1e-6.
- **Shorter grids, lower frequencies, another ratio.** The full table reproduces WRT to 1e-4 on grids of 25 frequencies
  from 0.05 Hz, of 26 to 34 frequencies from 0.03453 Hz, and of 44 frequencies from 0.025 Hz (ratio 1.1), and of 44 and 50
  frequencies with a ratio of 1.07, all with 24 directions. On the grid of 25 frequencies, the whole chain was also run,
  with one seed of practice spectra.

## Limits

- **Deep water.** The table holds the deep-water loci, which are the same at every frequency only on a grid with a constant
  frequency ratio. In finite depth, NL6 applies the WAM depth scaling, as WAVEWATCH III does for DIA and for WRT with
  `IQTYPE = 2`.
- **Grid length.** The loci of WRT reach about 12 to 14 bins below k1 and 20 bins above at a ratio of 1.1, and further
  with a smaller ratio or more directions. `build_table.py` takes their geometry from a bin that holds them, and it stops
  with a message if the grid is too short. Grids of 20 frequencies at a ratio of 1.1 and of 28 at 1.07 were too short.
- **Constants at low frequencies.** WRT computes its loci at its deep-water depth of 2000 m. The paper took k, cg and dk at
  1000 m, which differs from them by 2e-4 at 0.0345 Hz and more at lower frequencies. `make_operator.py` therefore uses
  Xnl's own values for grids that start below 0.033 Hz (`--constants auto`). `--constants paper` or `xnl` chooses either.
- **Directions over the full circle.** Xnl is run with `iq_grid = 3`. Directional sectors are not supported.
- **WRT's own settings.** The operator approximates Xnl with its defaults: 30 points per locus, no lumping, no
  Gauss–Legendre integration, and the run-time window of `k3/k1 <= 2.5`, `|theta3 - theta1| <= 75 deg`, and the density
  filter. `xnl_ref.py` checks these settings and stops if an `xnl4v5.cfg` would change them. The window has to be the same
  at every frequency and direction, which `build_table.py` checks.
- **Gravity.** Xnl is called with g = 9.81 m s-2, as in the paper. WAVEWATCH III calls its `NL2` with its own `GRAV` of
  9.806 m s-2, which changes the weights of the database by 0.25 %.
- **Build of Xnl.** `build_xnl.sh` compiles with `-O3`, as WAVEWATCH III's release build does, and ignores an `FFLAGS` of
  the environment (`XNL_FFLAGS` replaces `-O3`). The number of points that Xnl keeps on a few loci depends on rounding, and
  with `-O0` or `-O2` four of the 10 660 loci of the paper's grid differ. Another compiler may do the same. The toolchain is
  written to `tools/build/toolchain.txt`, and for the paper's grid `make_operator.py` warns if the database differs from
  the paper's.
- **Practice spectra.** Their peak frequencies (0.05 to 0.36 Hz) suit grids that reach from about 0.035 Hz to 1 Hz or
  higher. For another range, set `--fp-range` and `--swell-range` so that the peaks lie well inside the grid.
  `make_synthetic.py` warns if they do not.
