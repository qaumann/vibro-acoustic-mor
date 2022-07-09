% Produce the plots from the manuscript and write data to csv

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

%% set paths
dir_raw = 'output';
dir_graphics = [dir_raw filesep 'graphics'];
dir_data = [dir_raw filesep 'data'];

if exist(dir_graphics,'dir') ~= 7; mkdir(dir_graphics); end
if exist(dir_data,'dir') ~= 7; mkdir(dir_data); end

%% Figure 2: MORscores plate with hysteretic damping
load([dir_raw filesep 'data_plate_48_hysteretic'])
proj_methods = {'osimaginput', 'osrealinput'};

morscores = {};
for ii = 1:length(proj_methods)
    ind = contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    morscores{ii} = [this_results.morscore_linf]'; %#ok<SAGROW>
end

bar_names = {results.name};
for ii = 1:length(proj_methods)
    bar_names = strrep(bar_names, ['_' proj_methods{ii}],'');
end
bar_names = unique(strrep(bar_names, '_', ' '));

bar_cats = categorical(unique(bar_names));
bar_cats = reordercats(bar_cats, {'strint equi', 'strint avg std', ...
    'strint avg soa', 'strint avg sp', 'strint linf std', ...
    'strint linf soa', 'strint linf sp', 'minrel std', ...
    'minrel soa', 'minrel sp'});

figure('name','Manuscript Figure 2')
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'plate_48_hysteretic_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure 3: Plate with hysteretic damping: Relative L_inf errors
proj_methods = {'osimaginput'};
methods = {'avg_std', 'avg_soa', 'linf_std', 'linf_soa', 'minrel_std'};
errors = {};

for ii=1:length(methods)
    ind = contains({results.name}, proj_methods) & contains({results.name}, methods{ii});
    this_results = results(ind);
    
    errors{ii} = [this_results.linfrelerr]'; %#ok<SAGROW>
end

figure('name','Manuscript Figure 3')
semilogy(results(1).r, cell2mat(errors))
legend(strrep(methods,'_',' '))

qcsv([dir_data filesep 'plate_48_hysteretic_errs.csv'], ...
    [results(1).r' cell2mat(errors)], ...
    'header', ['r,' sprintf('%s,',methods{1:end-1}) methods{end}]);

%% Figure 4: MORscores plate with proportional damping
load([dir_raw filesep 'data_plate_48_rayleigh'])
proj_methods = {'osimaginput', 'osrealinput'};

morscores = {};
for ii = 1:length(proj_methods)
    ind = contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    morscores{ii} = [this_results.morscore_linf]'; %#ok<SAGROW>
end

bar_names = {results.name};
for ii = 1:length(proj_methods)
    bar_names = strrep(bar_names, ['_' proj_methods{ii}],'');
end
bar_names = unique(strrep(bar_names, '_', ' '));

bar_cats = categorical(unique(bar_names));
bar_cats = reordercats(bar_cats, {'strint equi', 'strint avg std', ...
    'strint avg soa', 'strint avg sp', 'strint linf std', ...
    'strint linf soa', 'strint linf sp', 'minrel std', ...
    'minrel soa', 'minrel sp', 'sobt'});

figure('name','Manuscript Figure 4')
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'plate_48_rayleigh_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure 5: Plate with proportional damping: Relative L_inf errors
proj_methods = {'osimaginput'};
methods = {'avg_std', 'linf_std', 'minrel_std', 'avg_sp', 'linf_sp', 'minrel_sp', ...
    'avg_soa', 'linf_soa', 'minrel_soa', 'equi', 'sobt'};
errors = {};

for ii=1:length(methods)
    ind = contains({results.name}, proj_methods) & contains({results.name}, methods{ii});
    this_results = results(ind);
    
    errors{ii} = [this_results.linfrelerr]'; %#ok<SAGROW>
end

figure('name','Manuscript Figure 5')
semilogy(results(1).r, cell2mat(errors))
legend(strrep(methods,'_',' '))

qcsv([dir_data filesep 'plate_48_rayleigh_errs.csv'], ...
    [results(1).r' cell2mat(errors)], ...
    'header', ['r,' sprintf('%s,',methods{1:end-1}) methods{end}]);

%% Figure 6: Plate with proportional damping: Relative errors at some r
sys = load_model('plate_48_rayleigh',false);
freq = abs(sys.s)/2/pi;
proj_methods = {'osimaginput'};
methods = {'linf_soa'};
errors = {};
rs = [140 220 240 250];
ind = contains({results.name}, proj_methods) & contains({results.name}, methods);
this_results = results(ind);

for ii=1:length(rs)
    errors{ii} = [this_results.relerr{rs(ii)}]'; %#ok<SAGROW>
end

figure('name','Manuscript Figure 6')
semilogy(freq, cell2mat(errors))
legend(rs)

qcsv([dir_data filesep 'plate_48_rayleigh_tferrs.csv'], ...
    [freq' cell2mat(errors)], ...
    'header', ['freq,' sprintf('%d,',rs(1:end-1)) num2str(rs(end))]);

%% Figure 7: Plate with proportional damping (SISO): SOBT formulas
sys = load_model('plate_48_rayleigh_single',false);
freq = abs(sys.s)/2/pi;
res = abs(sys.res);
load([dir_raw filesep 'data_plate_48_rayleigh_single'])
proj_methods = {'v', 'fv', 'pv', 'vp', 'p', 'so'};
tfs = {};
errors = {};
rs = 250;

for ii=1:length(proj_methods)
    ind = strcmp({results.name}, ['sobt_' proj_methods{ii}]);
    this_results = results(ind);
    
    tfs{ii} = [this_results.res{rs}]'; %#ok<SAGROW>
    errors{ii} = [this_results.relerr{rs}]'; %#ok<SAGROW>
end

figure('name','Manuscript Figure 7')
subplot(2,1,1)
plot(freq, 10*log10([res' cell2mat(tfs)]./1e-9))
legend([{'reference'} proj_methods(:)'])
subplot(2,1,2)
semilogy(freq, cell2mat(errors))
legend(proj_methods)

qcsv([dir_data filesep 'plate_48_rayleigh_single_sobt_tferrs.csv'], ...
    [freq' 10*log10([res' cell2mat(tfs)]./1e-9) cell2mat(errors)], ...
    'header', ['freq,res,' sprintf('%s,',proj_methods{:}) ...
    sprintf('err_%s,', proj_methods{1:end-1}) 'err_' proj_methods{end}]);

%% Figure x8: MORscores plate with proportional damping (SISO)
proj_methods = {'tsimag', 'tsreal', 'osimaginput', 'osrealinput', ...
    'osimagoutput', 'osrealoutput'};

morscores = {};
for ii = 1:length(proj_methods)
    ind = contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    morscores{ii} = [this_results.morscore_linf]'; %#ok<SAGROW>
end

bar_names = {results.name};
for ii = 1:length(proj_methods)
    bar_names = strrep(bar_names, ['_' proj_methods{ii}],'');
end
bar_names = unique(strrep(bar_names, '_', ' '));
bar_names(strcmp(bar_names,'sobt fv')) = [];
bar_names(strcmp(bar_names,'sobt p')) = [];
bar_names(strcmp(bar_names,'sobt pv')) = [];
bar_names(strcmp(bar_names,'sobt so')) = [];
bar_names(strcmp(bar_names,'sobt v')) = [];
bar_names(strcmp(bar_names,'sobt vp')) = [];

bar_cats = categorical(unique(bar_names));
bar_cats = reordercats(bar_cats, {'strint equi', 'strint avg std', ...
    'strint avg soa', 'strint avg sp', 'strint linf std', ...
    'strint linf soa', 'strint linf sp', 'minrel std', ...
    'minrel soa', 'minrel sp'});

figure('name','Manuscript Figure x8')
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'plate_48_rayleigh_single_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure x10: MORscores transmission
load([dir_raw filesep 'data_transmission'])
proj_methods = {'tsimag', 'tsreal', 'osimaginput', 'osrealinput', ...
    'osimagoutput', 'osrealoutput'};

morscores = {};
for ii = 1:length(proj_methods)
    ind = contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    morscores{ii} = [this_results.morscore_linf]'; %#ok<SAGROW>
end

bar_names = {results.name};
for ii = 1:length(proj_methods)
    bar_names = strrep(bar_names, ['_' proj_methods{ii}],'');
end
bar_names = unique(strrep(bar_names, '_', ' '));

bar_cats = categorical(unique(bar_names));
bar_cats = reordercats(bar_cats, {'strint equi', 'strint avg std', ...
    'strint avg soa', 'strint avg sp', 'strint linf std', ...
    'strint linf soa', 'strint linf sp', 'minrel std', ...
    'minrel soa', 'minrel sp'});

figure('name','Manuscript Figure x10')
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.6])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'transmission_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure x11: Plate with hysteretic damping: Relative L_inf errors
proj_methods = {'tsimag', 'tsreal', 'osimaginput', 'osrealinput', ...
    'osimagoutput', 'osrealoutput'};
methods = {'strint_equi'};
errors = {};

for ii=1:length(proj_methods)
    ind = contains({results.name}, methods) & contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    errors{ii} = [this_results.linfrelerr]'; %#ok<SAGROW>
end

figure('name','Manuscript Figure x11')
semilogy(this_results(1).r, cell2mat(errors))
legend(proj_methods)

qcsv([dir_data filesep 'transmission_errs.csv'], ...
    [this_results(1).r' cell2mat(errors)], ...
    'header', ['r,' sprintf('%s,',proj_methods{1:end-1}) proj_methods{end}]);

%% Figure x13: MORscores radiation
load([dir_raw filesep 'data_radiation'])
proj_methods = {'tsimag', 'tsreal', 'osimaginput', 'osrealinput', ...
    'osimagoutput', 'osrealoutput'};

morscores = {};
for ii = 1:length(proj_methods)
    ind = contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    morscores{ii} = [this_results.morscore_linf]'; %#ok<SAGROW>
end

bar_names = {results.name};
for ii = 1:length(proj_methods)
    bar_names = strrep(bar_names, ['_' proj_methods{ii}],'');
end
bar_names = unique(strrep(bar_names, '_', ' '));

bar_cats = categorical(unique(bar_names));
bar_cats = reordercats(bar_cats, {'strint equi', 'strint avg std', ...
    'strint avg soa', 'strint avg sp', 'strint linf std', ...
    'strint linf soa', 'strint linf sp', 'minrel std', ...
    'minrel soa', 'minrel sp'});

figure('name','Manuscript Figure x13')
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'radiation_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);













%% Helper function
function this_names = typeset_names(this_names)
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
