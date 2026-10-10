/* Flush subnormal floating-point numbers to zero (FTZ) and treat subnormal inputs as zero (DAZ) for the whole process: loaded with
   LD_PRELOAD=ftz/ftz.so by time_runs.py and count_steps.py --ftz (the source terms timed outside the model, bench/bench.py, run the same
   way). Build: gcc -O2 -shared -fPIC -o ftz/ftz.so ftz/ftz.c */
#include <xmmintrin.h>
#include <pmmintrin.h>
__attribute__((constructor)) static void set_ftz_daz(void) { _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON); _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON); }
