%% TEST_RUNME_PRESAMPLE
% Computes presampling data for the test model.

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

clear; close all

%% test model
sys = load_model('test', true);
wmin = 1*2*pi;
wmax = 1e4*2*pi;
ns = 50;
nsmin = 1;
nsmax = ns;

% standard presample
presample(sys, wmin, wmax, ns, 1, ns, []);

% multiple derivative presampling
ns = 17;
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns, []);

% second-order Arnoldi presample
r = 10;
ns = 5;
w_pre = 2*pi*[333 1024];
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns, w_pre);

