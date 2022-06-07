%% RUNME_PRESAMPLE
% Computes presampling data for all full order models.

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% model plate_48_rayleigh
sys = load_model('plate_48_rayleigh', true);
wmin = 1*2*pi;
wmax = 250*2*pi;
ns = 250;
side = 'input';

% standard presample
presample(sys, wmin, wmax, ns, 1, ns, [], side);

% multiple derivative presampling
ns = 84;
w_pre = 2*pi*[46 47 48 50];
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns, w_pre, side);

% second-order Arnoldi presample
r = 10;
ns = 25;
w_pre = 2*pi*[46 47 48 50];
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns, w_pre, side);

%% model plate_48_rayleigh_single
sys = load_model('plate_48_rayleigh_single', true);
wmin = 1*2*pi;
wmax = 250*2*pi;
ns = 250;

% standard presample
presample(sys, wmin, wmax, ns, 1, ns);

% multiple derivative presampling
ns = 84;
w_pre = 2*pi*[46 47 48 50];
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns, w_pre);

% second-order Arnoldi presample
r = 10;
ns = 25;
w_pre = 2*pi*[46 47 48 50];
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns, w_pre);

%% model plate_48_hysteretic
sys = load_model('plate_48_hysteretic', true);
wmin = 1*2*pi;
wmax = 250*2*pi;
ns = 250;
side = 'input';

% standard presample
presample(sys, wmin, wmax, ns, 1, ns, [], side);

% multiple derivative presampling
ns = 84;
w_pre = 2*pi*[46 47 48];
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns, w_pre, side);

% second-order Arnoldi presample
r = 10;
ns = 25;
w_pre = 2*pi*[46 47 48 50];
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns, w_pre, side);

%% model transmission
sys = load_model('transmission', true);
wmin = 1*2*pi;
wmax = 1000*2*pi;
ns = 200;

% standard presample
presample(sys, wmin, wmax, ns, 1, ns);

% multiple derivative presampling
ns = 67;
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns);

% second-order Arnoldi presample
r = 10;
ns = 20;
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns);

%% model radiation
sys = load_model('radiation', true);
wmin = 1*2*pi;
wmax = 600*2*pi;
ns = 200;

% standard presample
presample(sys, wmin, wmax, ns, 1, ns);

% multiple derivative presampling
ns = 67;
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns);

% second-order Arnoldi presample
r = 5;
ns = 40;
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns);

%% model poroacoustic
sys = load_model('poroacoustic', true);
wmin = 100*2*pi;
wmax = 1000*2*pi;
ns = 200;

% standard presample
presample(sys, wmin, wmax, ns, 1, ns);

% multiple derivative presampling
ns = 29;
presample_structure_preserving(sys, wmin, wmax, ns, 1, ns);

% second-order Arnoldi presample
r = 10;
ns = 20;
presample_aaa_arnoldi(sys, r, wmin, wmax, ns, 1, ns);
