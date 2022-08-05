Code, Data and Results for Numerical Experiments in "Structured model order
reduction for vibro-acoustic problems using interpolation and balancing
methods"
===========================================================================

This archive contains the companion codes and computed results for
the paper:

    Q. Aumann, S. W. R. Werner; Structured model order reduction for vibro-
    acoustic problems using interpolation and balancing methods,

which implement the reported numerical experiments.


## CODE

### Dependencies and installation

The computations were originally performed in MATLAB 9.9.0.1467703 (R2020b)
on single nodes of the MPI Magdeburg's computing cluster Mechthild running on
CentOS Linux version 7.9.2009 equipped with two 8 core Intel(R) Xeon(R)
Silver 4110 (Skylake) CPUs with a maximum clock rate of 3.0 GHz and 192 GB main
memory available.

Additionally, the following software packages need to be downloaded:

* M-M.E.S.S. version 2.0.1, doi:[10.5281/zenodo.3606345](https://doi.org/10.5281/zenodo.3606345)
* SOLBT version 3.0, doi:[10.5281/zenodo.4600763](https://doi.org/10.5281/zenodo.4600763)
* TOAR.m, by Ding Lu, available at [http://www.unige.ch/~dlu/toar.html](http://www.unige.ch/~dlu/toar.html)

which should be placed into the `software` directory for installation.

To properly set the MATLAB search path for the experiments, run the
`install.m` script before the experiments.


### Content and structure

The archive contains the following subdirectories:

* `models`: contains the benchmark models.
* `output`: contains a file for each benchmark with the raw result data compiled
  to one file.
* `presampling`: contains the data resulting from presampling strategies.
* `results`: contains all the computed results of the numerical
  experiments. Re-running the experiments will overwrite the content.
* `slurmscripts`: contains SLURM scripts to run experiments on the
  CoolMUC-2 computing cluster of the Leibniz supercomputing centre.
* `software`: location for additional software packages such as
  [M-M.E.S.S.](https://doi.org/10.5281/zenodo.3606345), 
  [SOLBT](https://doi.org/10.5281/zenodo.4600763), and
  [TOAR.m](http://www.unige.ch/~dlu/toar.html).
* `subroutines`: contains implementations and subroutines for the model
  order reduction strategies presented in the manuscript.

To test the computational setup before running the computational intensive
experiments, the `TEST_RUNME_*` scripts can be used.
The order in which to run these scripts is the same as for the later
experiments, with:

1. `TEST_RUNME_presample.m`: to test the presampling setup.
2. `TEST_RUNME_compute.m`: to test the model reduction setup.
3. `TEST_RUNME_evaluate_results.m`: to generate plots and save the results
    into .mat files.

The scripts to run the experiments lie in the same main directory as this
README and begin with `RUNME_*`. They perform the following experiments:

* `RUNME_presample.m`: Run the presampling for all numerical experiments.
* `RUNME_compute_[bench].m`: Run the model order reduction algorithms,
  where `[bench]` stands for one of the following experiments: 
  `plate_48_hysteretic`, `plate_48_rayleigh`, `plate_48_rayleigh_single`,
  `poroacoustic`, `radiation`, `transmission`.
* `RUNME_preprocess_results.m`: Compiles the raw result data from the folder 
  `results` to single files per benchmark. Produces the content of the folder
  `output`.
* `RUNME_plot_data.m`: Produces the plots shown in the manuscript.
* `RUNME_generate_data_tables.m`: Produces the data for the tables presented in
  the manuscript.

The scripts to run the experiments on an HPC cluster utilizing a SLURM
queue lie in the subdirectory `slurmscripts`. They perform the following
experiments:

* `batch_presample_[bench].sh`: Run the presampling algorithms, where
  `[bench]` stands for one of the following experiments: 
  `plate_48_hysteretic`, `plate_48_rayleigh`, `plate_48_rayleigh_single`,
  `poroacoustic`, `radiation`, `transmission`.
* `batch_compute_[bench].sh`: Run the model order reduction algorithms,
  where `[bench]` stands for one of the following experiments: 
  `plate_48_hysteretic`, `plate_48_rayleigh`, `plate_48_rayleigh_single`,
  `poroacoustic`, `radiation`, `transmission`.
* The settings for the SLURM environment can be modified in the scripts 
  `compute.sh`, `presample.sh`, `presample_aaa_arnoldi.sh`,
  `presample_structure_preserving.sh`, which are called from the `batch`
  scripts.

The experiments should be run in the order as given above, starting with
the `presample`, followed by all `compute` scripts and `evalute_results`
at the end.


## AUTHORS

Quirin Aumann
- affiliation: Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg (Germany)
- email: aumann@mpi-magdeburg.mpg.de
- orcid: [0000-0001-7942-5703](https://orcid.org/0000-0001-7942-5703)

Steffen W. R. Werner
- affiliation: Courant Institute of Mathematical Sciences, New York University (USA)
- email: steffen.werner@nyu.edu
- orcid: [0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)


## LICENSE

Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner

The software is licensed under the BSD-2-Clause License.
See [COPYING](COPYING) for a copy of the license.

The model and result data in the folders `models`, `presampling` and `result` is
licensed under a Creative Commons Attribution 4.0 License (CC-BY 4.0).
See [COPYING_DATA](COPYING_DATA) for a copy of the license.


## CITATION

### DOI

The DOI for version 1.1 of this archive is
[10.5281/zenodo.6806016](https://doi.org/10.5281/zenodo.6806016)


### Cite as

Q. Aumann and S. W. R. Werner. Code, data and results for numerical
experiments in "Structured model order reduction for vibro-acoustic
problems using interpolation and balancing methods" (version 1.1),
January 2022. doi:10.5281/zenodo.6806016


### BibTeX

    @MISC{supAumW22,
      author = {Aumann, Q. and Werner, S.~W.~R.},
      title  = {Code, Data and Results for Numerical Experiments in
                ``{S}tructured model order reduction for vibro-acoustic
                problems using interpolation and balancing methods''
                (version 1.1)},
      month  = aug,
      year   = {2022},
      doi    = {10.5281/zenodo.6806016}
    }
