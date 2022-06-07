# Setup OpenMP 
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads

# load the modules system 
source /etc/profile.d/modules.sh

# load matlab
module load apps/matlab/2020b
export MATLABPATH=.
export MKL_ENABLE_INSTRUCTIONS=AVX2;
