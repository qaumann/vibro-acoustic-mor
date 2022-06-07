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

echo "Starting presampling ${premethod} benchmark ${bench}:"
echo "presample_structure_preserving(load_model('$bench', true), $wmin, $wmax, $ns, 1, $ns, $wpre, '$side');"
echo

# Run MATLAB
matlab -nodisplay -r "install;\
sys=load_model('$bench', true);\
presample_structure_preserving(sys, $wmin, $wmax, $ns, 1, $ns, $wpre, '$side'); \
quit();"  < /dev/null


