#!/bin/bash
# Run directory for the 14-day global run (1-15 Jan 2025) with SWRT36 (default) or another operator (first argument, e.g.
# inputs/swrt9.ww3op for SWRT9), WAM depth scaling, semi-implicit integration with a time step of 15-900 s.
# Run ww3_shel with WW3_SWRT_TOPM=2 WW3_SWRT_TOPM_UST=1 (the weak-wind rule of the paper), see README.md.
set -euo pipefail
cd "$(dirname "$0")/.."
PKG="$PWD"; OP="$(realpath "${1:-inputs/swrt_K36.ww3op}")"; T="$(basename "${OP%.ww3op}")"; RUN="$PKG/run_${T}_global_14day"
test -x "$PKG/build/bin/ww3_shel" || { echo "build first" >&2; exit 1; }
test ! -e "$RUN" || { echo "refusing to overwrite $RUN" >&2; exit 1; }
mkdir -p "$RUN"; cd "$RUN"
ln -s "$PKG/inputs/Global-real.depth" .
ln -s "$PKG/inputs/wind_1deg_global_20250101to15.nc" .
cp "$PKG/inputs/restart_ic_20250101.ww3" restart.ww3
cp "$OP" swrt_operator.ww3op
cp "$PKG/inputs/points.list" points.list
ln -s "$PKG/switch" switch
cp "$PKG"/namelists/*.nml .
echo "run dir ready: $RUN"; ls -la "$RUN"
