% Read result files and arrange them in a struct

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

dir_results = 'results';
dir_presampling = 'presampling';
dir_output = 'output';

%% collect data
for bb = 1:length(benchs)
    
    bench = benchs{bb};
    
    % parameters for MOR score
    if strcmp(bench, 'radiation')
        morscore_eps = 1e-16;
        morscore_max = 200;
    elseif strcmp(bench, 'transmission')
        morscore_eps = 1e-16;
        morscore_max = 100;
    elseif strcmp(bench, 'poroacoustic')
        morscore_eps = 1e-16;
        morscore_max = 100;
    elseif strcmp(bench, 'plate_48_rayleigh')
        morscore_eps = 1e-16;
        morscore_max = 250;
    elseif strcmp(bench, 'plate_48_hysteretic')
        morscore_eps = 1e-6;
        morscore_max = 250;
    elseif strcmp(bench, 'plate_48_rayleigh_single')
        morscore_eps = 1e-16;
        morscore_max = 250;
    end
    
    % load full system response
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
        %distinguish plate_48_rayleigh from plate_48_rayleigh_single
        if strcmp(bench, 'plate_48_rayleigh') && contains(files(ii).name, 'single')
            continue
        end
        name = erase(files(ii).name, [bench + "_", ".mat"]);
        if startsWith(name, 'gram'); continue; end
        if startsWith(name, 'presampling'); continue; end
        if ~(contains(name,'equi') || contains(name,'aaaa') || ...
                contains(name,'strprs') || contains(name,'sobt'))
            tmp = strsplit(name,'_');
            name = [sprintf('%s_',tmp{1:end-1}),'std_',tmp{end}];
        end
        if contains(name,'osinput')
            name = strrep(name,'osinput','osimaginput');
        end
        results(ind).name = name;
        
        % load result
        load([files(ii).folder '/' files(ii).name]);
        tmp_r = [result_data.r];
        results(ind).r = 1:length(result_data);
        
        % compute relative errors per order
        
        for r = tmp_r
            if xor(all(imag(sys.res) == 0), all(imag(result_data(r).res) == 0))
                relerr = abs(abs(sys.res) - ...
                    abs(result_data(r).res(1:length(sys.res)))) ...
                    ./ abs(sys.res);
            else
                relerr = abs(sys.res - result_data(r).res(1:length(sys.res))) ...
                    ./ abs(sys.res);
            end
            results(ind).relerr{r} = relerr;
        end
        
        
        results(ind).maxrelerr = zeros(1,length(result_data));
        results(ind).hinfrelerr = zeros(1,length(result_data));
        
        % compute error norms. If full order result is abs, compute abs for ROM
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
        
        % fill error norms for not computed reduced orders
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
        
        % MOR computation time
        results(ind).ctime_mor = result_data(end).ctime.mor;
        
        % presampling computation time
        if contains(name,'aaaa')
            tmp = load([dir_presampling filesep 'presampling_aaaa_' bench '.mat'],'ctime');
            results(ind).ctime_presampling = tmp.ctime;
        elseif contains(name,'strprs')
            tmp = load([dir_presampling filesep 'presampling_strprs_' bench '.mat'],'ctime');
            results(ind).ctime_presampling = tmp.ctime;
        elseif contains(name,'std')
            tmp = load([dir_presampling filesep 'presampling_' bench '.mat'],'ctime');
            results(ind).ctime_presampling = tmp.ctime;
        elseif contains(name,'sobt_osimaginput')
            tmp = load([dir_presampling filesep bench '_gram_c.mat'],'ctime_gram_c');
            results(ind).ctime_presampling = tmp.ctime_gram_c;
        elseif contains(name,'sobt')
            tmp = load([dir_presampling filesep bench '_gram_c.mat'],'ctime_gram_c');
            results(ind).ctime_presampling = tmp.ctime_gram_c;
            tmp = load([dir_presampling filesep bench '_gram_o.mat'],'ctime_gram_o');
            results(ind).ctime_presampling = results(ind).ctime_presampling + tmp.ctime_gram_o;
        else
            results(ind).ctime_presampling = 0;
        end
        
        ind = ind + 1;
    end
    
    % add dummy for single-sided bt (real-valued not applicable)
    if any(strcmp({results.name},'sobt_osimaginput'))
        results(end+1).name = 'sobt_osrealinput';  %#ok<SAGROW>
        results(end).ctime_mor = 0;
        results(end).morscore_hinf = 0;
        results(end).morscore_maxrel = 0;
    end
    
    % sort results
    [~,idx]=sort({results.name});
    results = results(idx);
    
    if exist(dir_output, 'dir') ~= 7; mkdir(dir_output); end
    save([dir_output filesep 'data_' bench], 'results')
    
end
