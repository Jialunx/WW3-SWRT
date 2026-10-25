#!/bin/bash
#SBATCH --job-name=swrtk36
#SBATCH --partition=short,medium,long
#SBATCH --qos=priority
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# Submit from the package directory: cd /path/to/WW3-SWRT && sbatch scripts/submit_htc.slurm
set -euo pipefail
PKG="${SLURM_SUBMIT_DIR}"; BIN="$PKG/build/bin"; RUN="$PKG/run_swrtk36_global_14day"
source /etc/profile.d/00-modulepath.sh; source /etc/profile.d/modules.sh
module purge; module load foss/2023b netCDF-Fortran/4.6.1-gompi-2023b
cd "$RUN"
export WW3_SWRT_OPERATOR="$RUN/swrt_operator.ww3op"
echo "=== SWRT K=36 + WAM depth scaling, global 14 d, ARC HTC fast build, start $(date -u +%FT%TZ) on $(hostname), ntasks=${SLURM_NTASKS} ==="
grep -m1 "model name" /proc/cpuinfo
time "$BIN/ww3_grid" > ww3_grid.log 2>&1; echo "[ww3_grid done]"
time "$BIN/ww3_prnc" > ww3_prnc.log 2>&1; echo "[ww3_prnc done]"
echo "[ww3_shel starting]"
time srun --mpi=pmix -n "${SLURM_NTASKS}" --cpu-bind=cores --export=ALL,LD_PRELOAD="$PKG/build/ftz.so" "$BIN/ww3_shel" > ww3_shel.log 2>&1
echo "[ww3_shel done]"
time "$BIN/ww3_ounf" > ww3_ounf.log 2>&1; echo "[ww3_ounf done]"
cp ww3_ounp_spec.nml ww3_ounp.nml; time "$BIN/ww3_ounp" > ww3_ounp_spec.log 2>&1; echo "[ww3_ounp spectra done]"
cp ww3_ounp_src.nml ww3_ounp.nml; time "$BIN/ww3_ounp" > ww3_ounp_src.log 2>&1; echo "[ww3_ounp source terms done]"
grep -E "Initialization time|Elapsed time" ww3_shel.log | tee elapsed.txt
echo "=== finished at $(date -u +%FT%TZ) ==="
