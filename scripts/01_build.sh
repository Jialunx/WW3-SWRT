#!/bin/bash
# Build the depth-scaled sWRT (NL6) WAVEWATCH III on ARC HTC, fast version (2026-10-09):
# Release -O3 with AVX-512 for every CPU node of short/medium/long (Cascade Lake, Ice Lake login, Genoa),
# and ftz.so (flush subnormals to zero, as in the local timing runs). Run on htc-login.
set -euo pipefail
cd "$(dirname "$0")/.."
PKG="$PWD"; BUILD="$PKG/build"
source /etc/profile.d/00-modulepath.sh 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
module purge
module load foss/2023b netCDF-Fortran/4.6.1-gompi-2023b
module load CMake/3.27.6-GCCcore-13.2.0
module list 2>&1 | sed 's/^/    /'
rm -rf "$BUILD"; mkdir -p "$BUILD"; cd "$BUILD"
FL="-O3 -march=skylake-avx512 -mtune=cascadelake -DNDEBUG"
cmake "$PKG/ww3_src" -DSWITCH="$PKG/switch" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_FLAGS_RELEASE="$FL" -DCMAKE_C_FLAGS_RELEASE="$FL" 2>&1 | tee cmake.log
make -j"${MAKE_J:-8}" ww3_shel ww3_grid ww3_prnc ww3_ounf ww3_ounp 2>&1 | tee make.log
gcc -O2 -shared -fPIC -o "$BUILD/ftz.so" "$PKG/inputs/ftz.c"
grep -h "Fortran_FLAGS" "$BUILD/model/src/CMakeFiles/ww3_lib.dir/flags.make"
ls -la "$BUILD/bin" "$BUILD/ftz.so"
