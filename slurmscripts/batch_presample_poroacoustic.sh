#!/bin/bash

# BATCH_PRESAMPLE_POROACOUSTIC
# Queues presample scripts of the benchmark 'poroacoustic' using all applicable methods.

#
# This file is part of the Code, Data and Results for Numerical Experiments
# in "Structured model order reduction for vibro-acoustic problems using
# interpolation and balancing methods"
# Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
# All rights reserved.
# License: BSD 2-Clause License (see COPYING)
#

email="aumann@mpi-magdeburg.mpg.de"

bench='poroacoustic'
wmin=100*2*pi
wmax=1000*2*pi
ns=200
ns_aaaa=20
ns_strprs=29
r=10
side='ts'

sbatch -J pre_poro \
    --mail-user=${email} \
    --partition=long \
    --time=3-00:00:00 \
    --export=ALL,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns,wpre=[],side=$side \
    presample.sh

sbatch -J pre_poro_strprs \
    --mail-user=${email} \
    --partition=medium \
    --time=0-24:00:00 \
    --export=ALL,bench=$bench,wmin=$wmin,wmax=$wmax,ns=$ns_strprs,wpre=[],side=$side \
    presample_structure_preserving.sh

sbatch -J pre_poro_aaa \
    --mail-user=${email} \
    --partition=medium \
    --time=0-24:00:00 \
    --export=ALL,bench=$bench,r=$r,wmin=$wmin,wmax=$wmax,ns=$ns_aaaa,wpre=[],side=$side \
    presample_aaa_arnoldi.sh
