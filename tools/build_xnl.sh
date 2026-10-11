#!/bin/bash
# Build Van Vledder's Xnl, the WRT method of WAVEWATCH III's NL2, as the Python module xnl_wrt from the four Xnl files of
# ww3_src/model/src. The SWRT tools use it to make WRT's quadrature database for a grid, to read WRT's run-time constants,
# and to check the results against WRT. The files are compiled unchanged, except that two degree signs (Latin-1) in comments
# of mod_xnl4v5.f90 become '?' in the copies made here, since f2py reads its input as ASCII. mod_fileio.f90 takes only
# FILE_ENDIAN from WAVEWATCH III's CONSTANTS module; a one-line stand-in sets it to 'native' (bqf_decode.py reads either order).
# The flags are -O3, the optimization of WAVEWATCH III's release build: the number of points that Xnl keeps on a few loci
# depends on rounding (with -O0 or -O2, 4 of the 10660 loci of the paper's grid differ), and -O3 reproduces the database of
# the paper and of WAVEWATCH III. An FFLAGS of the environment is ignored; XNL_FFLAGS replaces -O3 if it is set.
# The compiler, NumPy and the flags are written to build/toolchain.txt.
# Needs gfortran and numpy (f2py); with Python >= 3.12, f2py also needs meson and ninja (pip install meson ninja).
# usage: bash tools/build_xnl.sh            ->  tools/build/xnl_wrt*.so, tools/build/toolchain.txt
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../ww3_src/model/src"
OUT="$HERE/build"
mkdir -p "$OUT/src"
cd "$OUT"
cat > src/constants_stub.f90 << 'EOF'
module constants
  implicit none
  character(len=20) :: file_endian = "native"
end module constants
EOF
FILES="src/constants_stub.f90"
for f in mod_constants mod_fileio serv_xnl4v5 mod_xnl4v5; do
  LC_ALL=C tr '\200-\377' '?' < "$SRC/$f.f90" > "src/$f.f90"
  FILES="$FILES src/$f.f90"
done
export FFLAGS="${XNL_FFLAGS:--O3}"
if python3 -m numpy.f2py -c $FILES -m xnl_wrt > build_xnl.log 2>&1; then
  { gfortran --version | head -1; python3 -c 'import sys, numpy; print("Python", sys.version.split()[0], "NumPy", numpy.__version__)'
    echo "FFLAGS=$FFLAGS"; } > toolchain.txt
  ls "$OUT"/xnl_wrt*.so
else
  echo "f2py failed, see $OUT/build_xnl.log"
  exit 1
fi
