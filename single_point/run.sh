#!/bin/bash
# Run one single-point simulation of the paper (see README.md in this directory).
#   bash single_point/run.sh <scenario> <bin> [operator]
#     <scenario>  g20, turn90, decay20, swell_cross, swell_opp or swell_wide
#     <bin>       bin/ directory of a single-point build: NL6 (SWRT), NL2 (WRT) or NL1 (DIA) in the switch file
#     [operator]  SWRT operator for an NL6 build (default: inputs/swrt_K30.ww3op)
# The run directory is single_point/run_<scenario>_<build directory>[_<operator>]. ww3_ounp writes the spectrum and
# the source terms every 30 min to ww3.202501_src.nc there.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
S="$1"; B="$(cd "$2" && pwd)"
[ -d "$HERE/$S" ] || { echo "unknown scenario: $S" >&2; exit 1; }
TAG="$(basename "$(dirname "$B")")"
if [ $# -ge 3 ]; then OP="$3"; TAG="${TAG}_$(basename "${OP%.ww3op}")"; else OP="$HERE/../inputs/swrt_K30.ww3op"; fi
OP="$(cd "$(dirname "$OP")" && pwd)/$(basename "$OP")"
D="$HERE/run_${S}_${TAG}"
rm -rf "$D"; mkdir -p "$D"; cp "$HERE/$S"/* "$D"/
cd "$D"
export WW3_SWRT_OPERATOR="$OP" OMP_NUM_THREADS=1
"$B/ww3_grid" > grid.log
"$B/ww3_strt" > strt.log
"$B/ww3_shel" > shel.log
"$B/ww3_ounp" > ounp.log
echo "done: $D/ww3.202501_src.nc"
