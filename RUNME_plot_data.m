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

fig = figure('name','Manuscript Figure 2');
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])
saveas(fig, [dir_graphics filesep 'figure_02.fig'])
saveas(fig, [dir_graphics filesep 'figure_02.png'])

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

fig = figure('name','Manuscript Figure 3');
semilogy(results(1).r, cell2mat(errors))
legend(strrep(methods,'_',' '))
saveas(fig, [dir_graphics filesep 'figure_03.fig'])
saveas(fig, [dir_graphics filesep 'figure_03.png'])

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

fig = figure('name','Manuscript Figure 4');
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])
saveas(fig, [dir_graphics filesep 'figure_04.fig'])
saveas(fig, [dir_graphics filesep 'figure_04.png'])

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

fig = figure('name','Manuscript Figure 5');
semilogy(results(1).r, cell2mat(errors))
legend(strrep(methods,'_',' '))
saveas(fig, [dir_graphics filesep 'figure_05.fig'])
saveas(fig, [dir_graphics filesep 'figure_05.png'])

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

fig = figure('name','Manuscript Figure 6');
semilogy(freq, cell2mat(errors))
legend(rs)
saveas(fig, [dir_graphics filesep 'figure_06.fig'])
saveas(fig, [dir_graphics filesep 'figure_06.png'])

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

fig = figure('name','Manuscript Figure 7');
subplot(2,1,1)
plot(freq, 10*log10([res' cell2mat(tfs)]./1e-9))
legend([{'reference'} proj_methods(:)'])
subplot(2,1,2)
semilogy(freq, cell2mat(errors))
legend(proj_methods)
saveas(fig, [dir_graphics filesep 'figure_02.fig'])
saveas(fig, [dir_graphics filesep 'figure_02.png'])

qcsv([dir_data filesep 'plate_48_rayleigh_single_sobt_tferrs.csv'], ...
    [freq' 10*log10([res' cell2mat(tfs)]./1e-9) cell2mat(errors)], ...
    'header', ['freq,res,' sprintf('%s,',proj_methods{:}) ...
    sprintf('err_%s,', proj_methods{1:end-1}) 'err_' proj_methods{end}]);

%% Figure 8: MORscores plate with proportional damping (SISO)
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

fig = figure('name','Manuscript Figure 8');
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])
saveas(fig, [dir_graphics filesep 'figure_08.fig'])
saveas(fig, [dir_graphics filesep 'figure_08.png'])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'plate_48_rayleigh_single_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure 10: MORscores transmission
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

fig = figure('name','Manuscript Figure x10');
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.6])
saveas(fig, [dir_graphics filesep 'figure_10.fig'])
saveas(fig, [dir_graphics filesep 'figure_10.png'])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'transmission_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure 11: Transmission: Relative L_inf errors
proj_methods = {'tsimag', 'tsreal', 'osimaginput', 'osrealinput', ...
    'osimagoutput', 'osrealoutput'};
methods = {'strint_equi'};
errors = {};

for ii=1:length(proj_methods)
    ind = contains({results.name}, methods) & contains({results.name}, proj_methods{ii});
    this_results = results(ind);
    
    errors{ii} = [this_results.linfrelerr]'; %#ok<SAGROW>
end

fig = figure('name','Manuscript Figure x11');
semilogy(this_results(1).r, cell2mat(errors))
legend(proj_methods)
saveas(fig, [dir_graphics filesep 'figure_11.fig'])
saveas(fig, [dir_graphics filesep 'figure_11.png'])

qcsv([dir_data filesep 'transmission_errs.csv'], ...
    [this_results(1).r' cell2mat(errors)], ...
    'header', ['r,' sprintf('%s,',proj_methods{1:end-1}) proj_methods{end}]);

%% Figure 13: MORscores radiation
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

fig = figure('name','Manuscript Figure 13');
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])
saveas(fig, [dir_graphics filesep 'figure_13.fig'])
saveas(fig, [dir_graphics filesep 'figure_13.png'])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'radiation_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure 14: Radiation: Relative L_inf errors
experiments = {'equi_tsimag', 'avg_std_tsimag', 'avg_soa_tsimag', ...
    'avg_soa_osimaginput', 'avg_soa_osimagoutput'}; 
errors = {};

for ii=1:length(experiments)
    ind = contains({results.name}, experiments{ii});
    this_results = results(ind);
    
    errors{ii} = [this_results.linfrelerr]'; %#ok<SAGROW>
end

fig = figure('name','Manuscript Figure 14');
semilogy(this_results(1).r, cell2mat(errors))
legend(strrep(experiments,'_',' '))
saveas(fig, [dir_graphics filesep 'figure_14.fig'])
saveas(fig, [dir_graphics filesep 'figure_14.png'])

qcsv([dir_data filesep 'radiation_errs.csv'], ...
    [this_results(1).r' cell2mat(errors)], ...
    'header', ['r,' sprintf('%s,',experiments{1:end-1}) experiments{end}]);

%% Figure 15: Radiation transfer functions
sys = load_model('radiation',false);
freq = abs(sys.s)/2/pi;
res = abs(sys.res(5,:));
experiments = {'equi_tsimag 140', 'equi_tsimag 200', 'avg_soa_tsimag 200', ...
    'avg_soa_osimagoutput 200'}; 
tfs = {};
errors = {};

for ii=1:length(experiments)
    tmp = strsplit(experiments{ii}, ' ');
    ind = contains({results.name}, tmp{1});
    this_results = results(ind);
    
    tfs{ii} = [this_results.res{str2double(tmp{2})}]'; %#ok<SAGROW>
    errors{ii} = [this_results.relerr{str2double(tmp{2})}]'; %#ok<SAGROW>
end

fig = figure('name','Manuscript Figure 15');
subplot(2,1,1)
plot(freq, 10*log10([res' cell2mat(tfs)]./1e-9))
xlim([0 600])
legend([{'reference'} strrep(experiments,'_',' ')])
subplot(2,1,2)
semilogy(freq, cell2mat(errors))
xlim([0 600])
legend(strrep(experiments,'_',' '))
saveas(fig, [dir_graphics filesep 'figure_15.fig'])
saveas(fig, [dir_graphics filesep 'figure_15.png'])

experiments = strrep(experiments,' ','_');
qcsv([dir_data filesep 'radiation_tferrs.csv'], ...
    [freq' 10*log10([res' cell2mat(tfs)]./1e-9) cell2mat(errors)], ...
    'header', ['freq,res,' sprintf('%s,',experiments{:}) ...
    sprintf('err_%s,', experiments{1:end-1}) 'err_' experiments{end}]);

%% Figure 17: MORscores poroacoustic
load([dir_raw filesep 'data_poroacoustic'])
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

fig = figure('name','Manuscript Figure 17');
bar(bar_cats, cell2mat(morscores))
legend(proj_methods)
ylim([0 0.4])
saveas(fig, [dir_graphics filesep 'figure_17.fig'])
saveas(fig, [dir_graphics filesep 'figure_17.png'])

mydat = [{typeset_names(bar_names)'} morscores(:)'];
qcsv([dir_data filesep 'poroacoustic_morscores.csv'], mydat, ...
    'header', ['method' sprintf(',%s',proj_methods{:})]);

%% Figure 18: Poroacoustic: Relative L_inf errors
proj_methods = {'tsimag'};
methods = {'avg_std', 'avg_sp', 'avg_soa', 'linf_std', 'linf_sp', 'linf_soa'};
errors = {};

for ii=1:length(methods)
    ind = contains({results.name}, methods{ii}) & contains({results.name}, proj_methods);
    this_results = results(ind);
    
    errors{ii} = [this_results.linfrelerr]'; %#ok<SAGROW>
end

figure('name','Manuscript Figure 18')
semilogy(this_results(1).r, cell2mat(errors))
legend(strrep(methods,'_',' '))
saveas(fig, [dir_graphics filesep 'figure_18.fig'])
saveas(fig, [dir_graphics filesep 'figure_18.png'])

qcsv([dir_data filesep 'poroacoustic_errs.csv'], ...
    [this_results(1).r' cell2mat(errors)], ...
    'header', ['r,' sprintf('%s,',methods{1:end-1}) methods{end}]);

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
