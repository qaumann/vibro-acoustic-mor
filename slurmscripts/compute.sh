#!/bin/bash
#SBATCH -o log/%j.%x.%N.out
#SBATCH -e log/%j.%x.%N.err
#SBATCH -D ..
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=RAM192
#SBATCH --partition=long
#SBATCH --mail-type=FAIL,END

source slurmscripts/load_modules.sh

echo "Starting $method $premethod benchmark ${bench}:"
echo "compute('$bench', '$method', '$pm', $ns, $rs, $wmin, $wmax, '$premethod', $wpre);"
echo

# Run MATLAB
matlab -nojvm -nosplash -batch "install;\
compute('$bench', '$method', '$pm', $ns, $rs, $wmin, $wmax, '$premethod', $wpre); \
quit();"  < /dev/null

