% Read preprocessed results and arrange them in tables

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

clear;
close all;

%% define parameters and output paths
benchs = {'plate_48_hysteretic', ...
    'plate_48_rayleigh', ...
    'plate_48_rayleigh_single', ...
    'transmission', ...
    'radiation', ...
    'poroacoustic'};

ordering = {'equi', 'avg_std', 'avg_sp', 'avg_soa', ...
    'linf_std', 'linf_sp', 'linf_soa', ...
    'minrel_std', 'minrel_sp', 'minrel_soa'};

dir_raw = 'output';
dir_data = [dir_raw filesep 'data'];

if exist(dir_data,'dir') ~= 7; mkdir(dir_data); end

%% collect data
for bb = 1:length(benchs)
    
    bench = benchs{bb};
    
    load([dir_raw filesep 'data_' bench])
    
    if strcmp(bench, 'plate_48_hysteretic') || strcmp(bench, 'plate_48_rayleigh')
        proj_methods = {'osimaginput', 'osrealinput'};
    else
        proj_methods = {'tsimag', 'tsreal', 'osimaginput', 'osrealinput', ...
            'osimagoutput', 'osrealoutput'};
    end
    
    if strcmp(bench, 'plate_48_rayleigh')
        ordering = {'equi', 'avg_std', 'avg_sp', 'avg_soa', ...
            'linf_std', 'linf_sp', 'linf_soa', ...
            'minrel_std', 'minrel_sp', 'minrel_soa', 'sobt'};
    else
        ordering = {'equi', 'avg_std', 'avg_sp', 'avg_soa', ...
            'linf_std', 'linf_sp', 'linf_soa', ...
            'minrel_std', 'minrel_sp', 'minrel_soa'};
    end
    
    ctimes = [];
    idx = 1;
    for pm = proj_methods
        for o = ordering
            ind = contains({results.name}, pm{:}) & contains({results.name}, o{:});
            ctimes(idx) = results(ind).ctime_mor;
            idx = idx + 1;
        end
    end
    
    ctimes = reshape(ctimes,length(ordering),length(proj_methods));
    T = splitvars(table(typeset_names(ordering'),ctimes));
    T = renamevars(T,1:width(T),[{'method'} proj_methods{:}]);
    
    qcsv([dir_data filesep 'table_' bench '_ctime.csv'], T);
    
    if strcmp(bench, 'plate_48_rayleigh_single')
        proj_methods = {'v', 'fv', 'pv', 'vp', 'p', 'so'};
        ctimes = zeros(1,length(proj_methods));
        idx = 1;
        for pm = proj_methods
            ind = strcmp({results.name}, ['sobt_' pm{:}]);
            ctimes(idx) = results(ind).ctime_mor;
            idx = idx+1;
        end
        T = splitvars(table({'\lmorbt{}'},ctimes));
        T = renamevars(T,1:width(T),[{'method'} proj_methods{:}]);
        
        qcsv([dir_data filesep 'table_' bench '_ctime_sobt.csv'], T);
    end
end

%% Helper function
function this_names = typeset_names(this_names)
this_names = strrep(this_names, '_', ' ');
this_names = strrep(this_names, 'strint ', '');
this_names = strrep(this_names, 'minrel', '\lmorminrel{}');
this_names = strrep(this_names, 'avg', '\lmoravg{}');
this_names = strrep(this_names, 'equi', '\lmorequi{}');
this_names = strrep(this_names, 'linf', '\lmorlinf{}');
this_names = strrep(this_names, 'sp', '\lmorsp{}');
this_names = strrep(this_names, 'soa', '\lmoraaa{}');
this_names = strrep(this_names, 'std', '\lmorstd{}');
this_names = strrep(this_names, 'sobt', '\lmorbt{}');
end
