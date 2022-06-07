% Evaluate test results and save them to csv and png

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

%% define parameters and paths
bench = 'test';

dir_results = 'results';
dir_graphics= 'results_graphics';
dir_data = 'results_data';

save_graphics_flag = true;
save_data_flag = true;

% parameters for MOR score
morscore_eps = 2.2204e-16;
morscore_max = 30;

%% load reduced order results and save them to a struct
sys = load_model(bench);
if size(sys.res,1) > 1
    sys.res = sys.res(5,:);
end

files = dir(dir_results);
results = struct();
ind = 1;

for ii=1:length(files)
    if ~startsWith(files(ii).name, bench); continue; end
    if ~endsWith(files(ii).name, '.mat'); continue; end
    name = erase(files(ii).name, [bench + "_", ".mat"]);
    if startsWith(name, 'gram'); continue; end
    if startsWith(name, 'presampling'); continue; end
    results(ind).name = name;
    
    % load result and compute error
    load([files(ii).folder '/' files(ii).name]);
    tmp_r = [result_data.r];
    results(ind).r = 1:length(result_data);
    results(ind).maxrelerr = zeros(1,length(result_data));
    results(ind).hinfrelerr = zeros(1,length(result_data));
    
    % compute errors. If full order result is abs, compute abs for ROM
    for jj=tmp_r
        this_result = result_data(jj).res;
        if length(sys.res) < length(this_result)
            this_result = this_result(1:length(sys.res));
        end
        
        if xor(all(imag(sys.res) == 0), all(imag(this_result) == 0))
            results(ind).maxrelerr(jj)  = max(abs(abs(sys.res) - abs(this_result)) ...
                ./ abs(sys.res));
            results(ind).hinfrelerr(jj) = max(abs(abs(sys.res) - abs(this_result))) ...
                ./ max(abs(sys.res));
        else
            results(ind).maxrelerr(jj)  = max(abs(sys.res - this_result) ...
                ./ abs(sys.res));
            results(ind).hinfrelerr(jj) = max(abs(sys.res - this_result)) ...
                ./ max(abs(sys.res));
        end
    end
    
    % fill errors for not computed reduced orders
    if length(tmp_r) ~= length(results(ind).r)
        prev_maxrelerr = 0;
        prev_hinfrelerr = 0;
        for jj=results(ind).r
            if isempty(result_data(jj).res)
                results(ind).maxrelerr(jj) = prev_maxrelerr;
                results(ind).hinfrelerr(jj) = prev_hinfrelerr;
            else
                prev_maxrelerr = results(ind).maxrelerr(jj);
                prev_hinfrelerr = results(ind).hinfrelerr(jj);
            end
        end
    end
    
    % fill missing errors with 1
    results(ind).maxrelerr(results(ind).maxrelerr == 0) = 1;
    results(ind).hinfrelerr(results(ind).hinfrelerr == 0) = 1;
    
    % truncate results for orders > morscore_max (from strprs and aaaa presampling)
    if length(results(ind).r) > morscore_max
        results(ind).r(morscore_max+1:end) = [];
        results(ind).maxrelerr(morscore_max+1:end) = [];
        results(ind).hinfrelerr(morscore_max+1:end) = [];
    end
    
    % compute MOR score
    tmp_hinfrelerr = [1 results(ind).hinfrelerr];
    tmp_hinfrelerr(tmp_hinfrelerr > 1) = 1;
    results(ind).morscore_hinf = trapz((0:morscore_max) / morscore_max, ...
        log10(tmp_hinfrelerr) / floor(log10(morscore_eps)));
    
    
    tmp_maxrelerr = [1 results(ind).maxrelerr];
    tmp_maxrelerr(tmp_maxrelerr > 1) = 1;
    results(ind).morscore_maxrel = trapz((0:morscore_max) / morscore_max, ...
        log10(tmp_maxrelerr) / floor(log10(morscore_eps)));
    
    ind = ind + 1;
end

% find available reduction methods
methods = {results.name};
methods = cellfun(@(x) erase(x, ...
    ["_tsimag", "_tsreal", "_osimaginput", "_osrealinput", ...
    "_osimagoutput", "_osrealoutput", "_osinput", "_osoutput", "_ts", ...
    "_fv", "_pv", "_pm", "_p", "_so", "_vpm", "_vp", "_v"]), ...
    methods, 'UniformOutput', false);
methods = unique(methods);
ind = find(strcmp(methods, 'structurereserving'));
if ~isempty(ind); methods{ind} = 'structure_preserving'; end

%% plot and export data
% figures: 
%           (1) max rel error over order per method
%           (2) hinf rel error over order per method
%           (3) max rel error over order for all methods (tsimag)
%           (4) hinf rel error over order for all methods (tsimag)
%           (5) tf and rel error for max(r) per method (tsimag)

if save_graphics_flag == true
    if exist(dir_graphics, 'dir') ~= 7; mkdir(dir_graphics);  end
end
if save_data_flag == true
    if exist(dir_data, 'dir') ~= 7; mkdir(dir_data);  end
end

%% max rel error over order per method
for m = methods
    fig = figure(1);
    set(fig,'defaulttextinterpreter','latex')
    
    selected_results = results(arrayfun(@(x) contains(x.name, m), results));
    maxr = max(arrayfun(@(x) max(x.r), selected_results));
    
    data = zeros(length(selected_results)+1, maxr);
    data(1,:) = 1:maxr;
    
    header = ['r,' strjoin({selected_results.name}, ',')];
    header = erase(header, [m{:} '_']);
    for ii=1:length(selected_results)
        name = erase(selected_results(ii).name, [m{:} '_']);
        semilogy(selected_results(ii).r, selected_results(ii).maxrelerr, ...
            'DisplayName', name)
        hold on
        data(ii+1,1:length(selected_results(ii).maxrelerr)) = selected_results(ii).maxrelerr;
    end
    hold off
    legend()
    title(['Pointwise error per order -- ' replace(m{:}, '_', ' ')])
    xlabel('Order $r$')
    ylabel('Maximum pointwise relative error')
    if save_graphics_flag == true
        saveas(gcf, [dir_graphics '/' bench '_maxrelerr_' m{:} '.png'], 'png');
        saveas(gcf, [dir_graphics '/' bench '_maxrelerr_' m{:} '.fig'], 'fig');
    end
    close(1)
    if save_data_flag == true
        qcsv([dir_data '/' bench '_maxrelerr_' m{:} '.csv'], transpose(data), ...
            'header', header);
    end
end

%% hinf rel error over order per method
for m = methods
    fig = figure(1);
    set(fig,'defaulttextinterpreter','latex')
    selected_results = results(arrayfun(@(x) contains(x.name, m), results));
    maxr = max(arrayfun(@(x) max(x.r), selected_results));
    
    data = zeros(length(selected_results)+1, maxr);
    data(1,:) = 1:maxr;
    
    header = ['r,' strjoin({selected_results.name}, ',')];
    header = erase(header, [m{:} '_']);
    for ii=1:length(selected_results)
        name = erase(selected_results(ii).name, [m{:} '_']);
        semilogy(selected_results(ii).r, selected_results(ii).hinfrelerr, ...
            'DisplayName', name)
        hold on
        data(ii+1,1:length(selected_results(ii).hinfrelerr)) = selected_results(ii).hinfrelerr;
    end
    hold off
    legend()
    title(['Relative $H_{\infty}$ error -- ' replace(m{:}, '_', ' ')])
    xlabel('Order $r$')
    ylabel('Relative $H_{\infty}$ error')
    if save_graphics_flag == true
        saveas(gcf, [dir_graphics '/' bench '_hinfrelerr_' m{:} '.png'], 'png');
        saveas(gcf, [dir_graphics '/' bench '_hinfrelerr_' m{:} '.fig'], 'fig');
    end
    close(1)
    if save_data_flag == true
        qcsv([dir_data '/' bench '_hinfrelerr_' m{:} '.csv'], transpose(data), ...
            'header', header);
    end
end

%% max rel error over order for all methods
if startsWith(sys.name, 'plate')
    proj_method = 'osimaginput';
else
    proj_method = 'tsimag';
end
selected_results = results(arrayfun(@(x) contains(x.name, proj_method), results));
maxr = max(arrayfun(@(x) max(x.r), selected_results));

data = zeros(length(selected_results)+1, maxr);
data(1,:) = 1:maxr;

header = ['r,' strjoin({selected_results.name}, ',')];
header = erase(header, ['_' proj_method]);
    
fig = figure(1);
set(fig,'defaulttextinterpreter','latex')
for ii=1:length(selected_results)
    name = erase(selected_results(ii).name, ['_' proj_method]);
    semilogy(selected_results(ii).r, selected_results(ii).maxrelerr, ...
        'DisplayName', replace(name, '_', ' '))
    hold on
    data(ii+1,1:length(selected_results(ii).maxrelerr)) = selected_results(ii).maxrelerr;
end
hold off
legend()
title('Pointwise error per order -- tsimag')
xlabel('Order $r$')
ylabel('Maximum pointwise relative error')
if save_graphics_flag == true
    saveas(gcf, [dir_graphics '/' bench '_maxrelerr_' proj_method '.png'], 'png');
    saveas(gcf, [dir_graphics '/' bench '_maxrelerr_' proj_method '.fig'], 'fig');
end
close(1)
if save_data_flag == true
    qcsv([dir_data '/' bench '_maxrelerr_' proj_method '.csv'], transpose(data), ...
        'header', header);
end


%% hinf rel error over order for all methods (tsimag)
if startsWith(sys.name, 'plate') && ~contains(sys.name, 'single')
    proj_method = 'osimaginput';
else
    proj_method = 'tsimag';
end
selected_results = results(arrayfun(@(x) contains(x.name, proj_method), results));
maxr = max(arrayfun(@(x) max(x.r), selected_results));

data = zeros(length(selected_results)+1, maxr);
data(1,:) = 1:maxr;

header = ['r,' strjoin({selected_results.name}, ',')];
header = erase(header, ['_' proj_method]);
    
fig = figure(1);
set(fig,'defaulttextinterpreter','latex')
for ii=1:length(selected_results)
    name = erase(selected_results(ii).name, ['_' proj_method]);
    semilogy(selected_results(ii).r, selected_results(ii).hinfrelerr, ...
        'DisplayName', replace(name, '_', ' '))
    hold on
    data(ii+1,1:length(selected_results(ii).hinfrelerr)) = selected_results(ii).hinfrelerr;
end
hold off
legend()
title(['Relative $H_{\infty}$ error -- ' proj_method])
xlabel('Order $r$')
ylabel('Relative $H_{\infty}$ error')
if save_graphics_flag == true
    saveas(gcf, [dir_graphics '/' bench '_hinfrelerr_' proj_method '.png'], 'png');
    saveas(gcf, [dir_graphics '/' bench '_hinfrelerr_' proj_method '.fig'], 'fig');
end
close(1)
if save_data_flag == true
    qcsv([dir_data '/' bench '_hinfrelerr_' proj_method '.csv'], transpose(data), ...
        'header', header);
end

%% tf and rel error for max(r) per method (tsimag)
for m = methods
    proj_methods = ["tsimag", "tsreal", "osimaginput", "osrealinput", ...
        "osimagoutput", "osrealoutput", "fv", "p", "pm", "pv", "so", "v", "vp", "vpm"];
    for pm = proj_methods
        if exist([dir_results '/' bench '_' m{:} '_' pm{:} '.mat'], 'file') == 2
            load([dir_results '/' bench '_' m{:} '_' pm{:} '.mat'], 'result_data');
        else 
            continue
        end
        fig = figure(1);
        set(fig,'defaulttextinterpreter','latex')
        subplot(2,1,1)
        semilogy(abs(sys.s)/2/pi, abs(sys.res), ...
            abs(sys.s)/2/pi, abs(result_data(end).res(1:length(sys.res))));
        xlabel('Frequency $\omega$ (rad/sec)');
        ylabel('Magnitude');
        legend({'Full model', ['$r=' num2str(result_data(morscore_max).r) '$']}, 'interpreter', 'latex');
        title(['Transfer functions -- ' replace(m{:}, '_', ' ')]);

        subplot(2, 1, 2);
        if isempty(result_data(morscore_max).res)
            ind = length(result_data);
        else
            ind = morscore_max;
        end

        if xor(all(imag(sys.res) == 0), all(imag(result_data(ind).res) == 0))
            hinferr = abs(abs(sys.res) - ...
                abs(result_data(ind).res(1:length(sys.res)))) ...
                ./ abs(sys.res);
            semilogy(abs(sys.s)/2/pi, hinferr);
        else
            hinferr = abs(sys.res - result_data(ind).res(1:length(sys.res))) ...
                ./ abs(sys.res);
            semilogy(abs(sys.s)/2/pi, hinferr);
        end

        xlabel('Frequency $\omega$ (rad/sec)');
        ylabel('Magnitude');
        title('Relative errors');

        if save_graphics_flag == true
            saveas(gcf, [dir_graphics '/' bench '_tf_' m{:} '_' pm{:} '.png'], 'png');
            saveas(gcf, [dir_graphics '/' bench '_tf_' m{:} '_' pm{:} '.fig'], 'fig');
        end
        close(1)

        if save_data_flag == true
            data = [abs(sys.s)/2/pi; abs(sys.res); ...
                abs(result_data(end).res(1:length(sys.res))); hinferr];
            qcsv([dir_data '/' bench '_tf_' m{:} '_' pm{:} '.csv'], transpose(data), ...
                'header', 's,res,res_r,hinferr' );
        end
    end
end

%% MOR score
% save MOR scores as csv with structure
% [method, presampling, projection, morscore_hinf, morscore_maxrel]
mor_methods = ["strint_avg", "strint_equi", "strint_linf", "minrel"];
pre_methods = ["aaaa", "strprs"];
proj_methods = ["tsimag", "tsreal", "osimaginput", "osrealinput", ...
    "osimagoutput", "osrealoutput"];

data_methods = cell(length(results), 1);
data_pre_methods = cell(length(results), 1);
data_proj_methods = cell(length(results), 1);
for ii=1:length(results)
    ind = arrayfun(@(x) contains(results(ii).name, x), mor_methods);
    data_methods{ii} = mor_methods(ind);
    ind = arrayfun(@(x) contains(results(ii).name, x), pre_methods);
    data_pre_methods{ii} = pre_methods(ind);
    ind = arrayfun(@(x) contains(results(ii).name, x), proj_methods);
    data_proj_methods{ii} = proj_methods(ind);
end
    
if save_data_flag == true
    qcsv([dir_data '/' bench '_mor_score.csv'], ...
        {data_methods, data_pre_methods, data_proj_methods, ...
        {results.morscore_hinf}', {results.morscore_maxrel}'}, ...
        'header', 'method,presampling,projection,morscore_hinf,morscore_maxrel');
end

% save MOR score per projection method
if save_data_flag == true
    proj_methods = ["tsimag", "tsreal", "osimaginput", "osrealinput", ...
        "osimagoutput", "osrealoutput"];
    
    for m = proj_methods
        ind = contains({results.name}, m{:});
        if ~any(ind); continue; end
        this_results = results(ind);
        
        % typeset names
        this_names = {this_results.name}';
        this_names = strrep(this_names, ['_' m{:}],'');
        this_names = strrep(this_names, '_', ' ');
        this_names = strrep(this_names, 'strint ', '');
        this_names = strrep(this_names, 'minrel', '\lmorminrel{}');
        this_names = strrep(this_names, 'avg', '\lmoravg{}');
        this_names = strrep(this_names, 'equi', '\lmorequi{}');
        this_names = strrep(this_names, 'linf', '\lmorlinf{}');
        this_names = strrep(this_names, 'strprs', '\lmorsp{}');
        this_names = strrep(this_names, 'aaaa', '\lmoraaa{}');
        this_names = strrep(this_names, 'pod', '\lmorpod{}');
        
        qcsv([dir_data '/' bench '_mor_score_' m{:} '.csv'], ...
            {this_names {this_results.morscore_hinf}' {this_results.morscore_maxrel}'}, ...
            'header', 'method,morscore_hinf,morscore_maxrel');
    end
end

