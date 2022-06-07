#!/bin/bash

# BATCH_PLATE_48_HYSTERETIC
# Queues compute scripts of the benchmark 'plate_48_hysteretic' using all applicable methods.

#
# This file is part of the Code, Data and Results for Numerical Experiments
# in "Structured model order reduction for vibro-acoustic problems using
# interpolation and balancing methods"
# Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
# All rights reserved.
# License: BSD 2-Clause License (see COPYING)
#

email="aumann@mpi-magdeburg.mpg.de"

bench='plate_48_hysteretic'
wmin=1*2*pi
wmax=250*2*pi
ns=250
ns_aaaa=25
ns_strprs=84
rs=1:250
rs_equi=5:250
rs_equi_real=10:250
w_pre_equi="2*pi*[44:48]"

for pm in 'osimaginput' 'osrealinput'; do
	for method in 'strint_avg' 'strint_linf' 'minrel'; do
		sbatch -J ph_${method}_${pm} \
            --mail-user=${email} \
            --partition=medium \
            --time=0-24:00:00 \
            --export=ALL,method=${method},bench=${bench},wmin=$wmin,wmax=$wmax,ns=$ns,rs=$rs,pm=$pm,wpre=[] \
            compute.sh

		sbatch -J ph_${method}_${pm} \
            --mail-user=${email} \
            --partition=medium \
            --time=0-24:00:00 \
            --export=ALL,method=${method},premethod=aaaa,bench=${bench},wmin=$wmin,wmax=$wmax,ns=$ns_aaaa,rs=$rs,pm=$pm,wpre=[] \
            compute.sh

		sbatch -J ph_${method}_${pm} \
            --mail-user=${email} \
            --partition=medium \
            --time=0-24:00:00 \
            --export=ALL,method=${method},premethod=strprs,bench=${bench},wmin=$wmin,wmax=$wmax,ns=$ns_strprs,rs=$rs,pm=$pm,wpre=[] \
            compute.sh
	done
done

method='strint_equi'
pm='osimaginput'
sbatch -J ph_${method}_${pm} \
    --mail-user=${email} \
    --partition=long \
    --time=7-00:00:00 \
    --export=ALL,method=${method},bench=${bench},wmin=$wmin,wmax=$wmax,ns=$ns,rs=$rs_equi,pm=$pm,wpre="$w_pre_equi" \
    compute.sh

pm='osrealinput'
sbatch -J ph_${method}_${pm} \
    --mail-user=${email} \
    --partition=long \
    --time=7-00:00:00 \
    --export=ALL,method=${method},bench=${bench},wmin=$wmin,wmax=$wmax,ns=$ns,rs=$rs_equi_real,pm=$pm,wpre="$w_pre_equi" \
    compute.sh
