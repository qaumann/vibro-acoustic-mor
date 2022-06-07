#!/bin/bash

# BATCH_PRESAMPLE_PLATE_48_HYSTERETIC
# Queues presample scripts of the benchmark 'plate_48_hysteretic' using all applicable methods.

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
r=10
w_pre_strprs="2*pi*[46:48]"
w_pre_aaaa="2*pi*[46 47 48 50]"
side='input'

sbatch -J pre_ph \
    --mail-user=${email} \
    --time=3-00:00:00 \
    --export=ALL,\
    bench=$bench,\
    wmin=$wmin,\
    wmax=$wmax,\
    ns=$ns,\
    wpre=[],\
    side=$side \
    presample.sh

sbatch -J pre_ph_strprs \
    --mail-user=${email} \
    --time=3-00:00:00 \
    --export=ALL,\
    bench=$bench,\
    wmin=$wmin,\
    wmax=$wmax,\
    ns=$ns_strprs,\
    wpre="$w_pre_strprs",\
    side=$side \
    presample_structure_preserving.sh

sbatch -J pre_ph_aaa \
    --mail-user=${email} \
    --time=3-00:00:00 \
    --export=ALL,\
    bench=$bench,\
    r=$r,\
    wmin=$wmin,\
    wmax=$wmax,\
    ns=$ns_aaaa,\
    wpre="$w_pre_aaaa",\
    side=$side \
    presample_aaa_arnoldi.sh
