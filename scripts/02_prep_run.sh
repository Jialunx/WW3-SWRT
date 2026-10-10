#!/bin/bash
# Run directory for the 14-day global run (1-15 Jan 2025), K = 30, WAM depth scaling; inputs identical to the local
# runs/swrtk30dscale_global_14day_ic20250101 (same namelists, restart_ic_20250101, ERA5 winds, operator).
set -euo pipefail
cd "$(dirname "$0")/.."
PKG="$PWD"; RUN="$PKG/run_swrtk30dscale_global_14day"
test -x "$PKG/build/bin/ww3_shel" || { echo "build first" >&2; exit 1; }
test ! -e "$RUN" || { echo "refusing to overwrite $RUN" >&2; exit 1; }
mkdir -p "$RUN"; cd "$RUN"
ln -s "$PKG/inputs/Global-real.depth" .
ln -s "$PKG/inputs/wind_1deg_global_20250101to15.nc" .
cp "$PKG/inputs/restart_ic_20250101.ww3" restart.ww3
cp "$PKG/inputs/swrt_K30.ww3op" swrt_operator.ww3op
cp "$PKG/inputs/points.list" points.list
ln -s "$PKG/switch" switch
cp "$PKG"/namelists/*.nml .
echo "run dir ready: $RUN"; ls -la "$RUN"
