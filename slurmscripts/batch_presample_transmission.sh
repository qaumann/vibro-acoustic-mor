#!/bin/bash

# BATCH_PRESAMPLE_TRANSMISSION
# Queues presample scripts of the benchmark 'transmission' using all applicable methods.

#
# This file is part of the Code, Data and Results for Numerical Experiments
# in "Structured model order reduction for vibro-acoustic problems using
# interpolation and balancing methods"
# Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
# All rights reserved.
# License: BSD 2-Clause License (see COPYING)
#

email="aumann@mpi-magdeburg.mpg.de"

bench='transmission'
wmin=1*2*pi
wmax=1000*2*pi
ns=200
ns_aaaa=20
ns_strprs=67
r=10
side='ts'

sbatch -J pre_tr \
    --mail-user=${email} \
    --partition=medium \
    --time=0-24:00:00 \
    --export=ALL,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns,wpre=[],side=$side \
    presample.sh

sbatch -J pre_tr_strprs \
    --mail-user=${email} \
    --partition=medium \
    --time=0-24:00:00 \
    --export=ALL,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns_strprs,wpre=[],side=$side \
    presample_structure_preserving.sh

sbatch -J pre_tr_aaa \
    --mail-user=${email} \
    --partition=medium \
    --time=0-24:00:00 \
    --export=ALL,bench=$bench,r=$r,wmin=$wmin,wmax=$wmax,ns=$ns_aaaa,wpre=[],side=$side \
    presample_aaa_arnoldi.sh
