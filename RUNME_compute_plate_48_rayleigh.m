%% RUNME_COMPUTE_PLATE_48_RAYLEIGH
% Computes reduced models of the benchmark 'plate_48_rayleigh' using all applicable methods.

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% settings
bench = 'plate_48_rayleigh';
wmin = 1*2*pi;
wmax = 250*2*pi;
ns = 250;
ns_aaaa = 25;
ns_strprs = 84;
rs = 1:250;
rs_equi = 5:250;
w_pre_equi = 2*pi*(44:48);

%% perform interpolation algorithms
methods = {'strint_avg', 'strint_linf', 'minrel'};
proj_methods = {'osimaginput', 'osrealinput'};
pre_methods = {'aaaa', 'strprs'};

for proj_method = proj_methods
    for method = methods
        compute(bench, ...
            method{:}, proj_method{:}, ...
            ns, rs, wmin, wmax);
        
        compute(bench, ...
            method{:}, proj_method{:}, ...
            ns_aaaa, rs, wmin, wmax, ...
            'aaaa');
        
        compute(bench, ...
            method{:}, proj_method{:}, ...
            ns_strprs, rs, wmin, wmax, ...
            'strprs');
    end
    
    compute(bench, ...
        'strint_equi', proj_method{:}, ...
        ns, rs_equi, wmin, wmax, ...
        '', w_pre_equi);
    
end
