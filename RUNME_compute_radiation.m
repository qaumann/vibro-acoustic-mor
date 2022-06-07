%% RUNME_COMPUTE_RADIATION
% Computes reduced models of the benchmark 'radiation' using all applicable methods.

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% settings
bench = 'radiation';
wmin = 1*2*pi;
wmax = 600*2*pi;
ns = 200;
ns_aaaa = 40;
ns_strprs = 67;
rs = 1:200;

%% perform interpolation algorithms
methods = {'strint_avg', 'strint_linf', 'minrel'};
proj_methods = {'tsimag', 'tsreal', ...
    'osimaginput', 'osimagoutput', 'osrealinput', 'osrealoutput'};
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
        ns, rs, wmin, wmax);
    
end
