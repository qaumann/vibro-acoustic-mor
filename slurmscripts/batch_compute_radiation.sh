#!/bin/bash

# BATCH_RADIATION
# Queues compute scripts of the benchmark 'radiation' using all applicable methods.

#
# This file is part of the Code, Data and Results for Numerical Experiments
# in "Structured model order reduction for vibro-acoustic problems using
# interpolation and balancing methods"
# Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
# All rights reserved.
# License: BSD 2-Clause License (see COPYING)
#

email="aumann@mpi-magdeburg.mpg.de"

bench='radiation'
wmin=1*2*pi
wmax=600*2*pi
ns=200
ns_aaaa=40
ns_strprs=67
rs=1:200

for pm in 'tsimag' 'tsreal' 'osimaginput' 'osimagoutput' 'osrealinput' 'osrealoutput'; do
	for method in 'strint_avg' 'strint_linf' 'minrel'; do
		sbatch -J rad_${method}_${pm} \
            --mail-user=${email} \
            --partition=long \
            --time=3-00:00:00 \
            --export=ALL,method=$method,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns,rs=$rs,pm=$pm,wpre=[] \
            compute.sh

		sbatch -J rad_${method}_${pm}_aaa \
            --mail-user=${email} \
            --partition=long \
            --time=3-00:00:00 \
            --export=ALL,method=$method,premethod=aaaa,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns_aaaa,rs=$rs,pm=$pm,wpre=[] \
            compute.sh

		sbatch -J rad_${method}_${pm}_strprs \
            --mail-user=${email} \
            --partition=long \
            --time=3-00:00:00 \
            --export=ALL,method=$method,premethod=strprs,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns_strprs,rs=$rs,pm=$pm,wpre=[] \
            compute.sh
	done

	method='strint_equi'
	sbatch -J rad_${method}_${pm} \
        --mail-user=${email} \
        --partition=long \
        --time=7-00:00:00 \
        --export=ALL,method=$method,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns,rs=$rs,pm=$pm,wpre=[] \
        compute.sh
done
