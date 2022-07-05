function [Ar, br, cr, s, res] = compute(bench, method, proj_method, ns, rs, wmin, wmax, pre_method, w_pre, r0)
%COMPUTE compute reduced order models using different reduction methods
%
% SYNTAX:
%   [Ar, br, cr, s, res] = COMPUTE(bench, method, proj_method, ns, rs, ...
%           wmin, wmax, pre_method, w_pre)
%
% DESCRIPTION:
%
% INPUT:
%   bench - full order model name
%   method - MOR method
%              'strint_avg'  - structured interpolation with approximate 
%                              subspaces
%              'strint_equi' - structured interpolation with equidistant 
%                              points
%              'strint_linf' - structured interpolation with L-infinity
%                              points
%              'minrel'      - minimal realization algorithm
%              'sobt'        - balanced truncation
%   proj_method - type of projection space for interpolation, resp. BT
%                 formula
%              'tsimag'       - two-sided projection allowing imag. parts
%              'tsreal'       - two-sided projection only real parts
%              'osimaginput'  - one-sided projection imag. input parts
%              'osrealinput'  - one-sided projection real input parts
%              'osimagoutput' - one-sided projection imag. output parts
%              'osrealoutput' - one-sided projection real output parts
%              'p'            - position balancing
%              'pm'           - position balancing (diagonalized M)
%              'pv'           - position-velocity balancing
%              'vp'           - velocity-position balancing
%              'vpm'          - velocity-position balancing (diag. M)
%              'v'            - velocity balancing
%              'fv'           - free velocity balancing
%              'so'           - second-order balancing
%   ns - number of frequency points
%   rs - array of sizes for the reduced models
%   wmin, wmax - min and max expansion frequency
%   pre_method - type of presampling method (for strint_equi, strint_linf,
%                minrel); may be empty
%              'strprs' - structure preserving approach by Beattie/Gugercin
%              'aaaa'   - second order Arnoldi with arbitrary interpolation
%                         order
%   w_pre - predefined frequency points; may be empty
%   r0 - Subspace order of Arnoldi presampling basis; may be empty

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP SYSTEM DATA.                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'Computing reduced model for benchmark %s:\n', bench);
if nargin == 7
    fprintf(1, '\tMethod: %s_%s\n', method, proj_method);
    pre_method = '';
elseif nargin == 8
    fprintf(1, '\tMethod: %s_%s_%s\n', method, pre_method, proj_method);
end

if nargin < 9
    w_pre = [];
end

fprintf('Compute process PID=%d running on node %s with %d threads.\n',...
    feature('getpid'),...
    strtrim(evalc('system(''hostname'');')),...
    maxNumCompThreads);


fprintf(1, 'Load sytem data.\n');
fprintf(1, '----------------\n');
ctime_all = tic;

sys = load_model(bench);
result = struct();

fprintf(1, '\n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD PRESAMPLING DATA.                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, 'Load presampling data.\n');
fprintf(1, '----------------------\n');

if strcmpi(method, 'strint_avg') || ...
        strcmp(method, 'strint_linf') || ...
        strcmp(method, 'minrel')
    
    if isempty(pre_method)
        fname = ['presampling/presampling_' sys.name '.mat'];
    else
        if exist('r0','var') == 1
            fname = ['presampling/presampling_' pre_method '_' sys.name ...
                '_r_' num2str(r0) '.mat'];
        else
            fname = ['presampling/presampling_' pre_method '_' sys.name '.mat'];
        end
    end
    
    if exist(fname, 'file') == 2
        tmp = load(fname);

        if ns ~= length(tmp.w)
            error('Number of interpolation points does not match')
        end

        if 1i*wmin ~= min(tmp.w) || 1i*wmax ~= max(tmp.w)
            error('Bounds for frequency range do not match')
        end
        w = tmp.w;

        if strcmpi(method, 'strint_avg') || strcmpi(method, 'minrel')
            if strcmpi(proj_method, 'tsimag') ...
                    || strcmpi(proj_method, 'osimaginput')
                V = tmp.V;
            end
            if strcmpi(proj_method, 'tsreal') ...
                    || strcmpi(proj_method, 'osrealinput')
                V = zeros(size(tmp.V,1), 2*size(tmp.V,2));
                V(:,1:2:end) = real(tmp.V);
                V(:,2:2:end) = imag(tmp.V);
            end
            if strcmpi(proj_method, 'tsimag') ...
                    || strcmpi(proj_method, 'osimagoutput')
                W = tmp.W;
            end
            if strcmpi(proj_method, 'tsreal') ...
                    || strcmpi(proj_method, 'osrealoutput')
                W = zeros(size(tmp.W,1), 2*size(tmp.W,2));
                W(:,1:2:end) = real(tmp.W);
                W(:,2:2:end) = imag(tmp.W);
            end
            
            if strcmpi(proj_method, 'osimaginput') ...
                    || strcmpi(proj_method, 'osrealinput')
                W = V;
            elseif strcmpi(proj_method, 'osimagoutput') ...
                    || strcmpi(proj_method, 'osrealoutput')
                V = W;
            end

            if strcmpi(method, 'minrel')
                % Compute dominant subspaces.
                time_minrel_subspaces = tic;

                lngthA = length(sys.A) - sum(cellfun(@isempty, sys.A));
                sizeV  = size(V, 2);

                if strcmpi(proj_method, 'tsimag') ...
                        || strcmpi(proj_method, 'tsreal') ...
                        || strcmpi(proj_method, 'osimaginput') ...
                        || strcmpi(proj_method, 'osrealinput')
                    WAVT = zeros(lngthA * sizeV, sizeV);

                    for k = 1:lngthA
                        WAVT((k-1)*sizeV+1:k*sizeV, :) = W' * (sys.A{k} * V);
                    end

                    [~, ~, Vtmp] = svd(WAVT, 'econ');

                    if strcmpi(proj_method, 'osimaginput') ...
                            || strcmpi(proj_method, 'osrealinput')
                        Wtmp = Vtmp;
                    end
                end

                if strcmpi(proj_method, 'tsimag') ...
                        || strcmpi(proj_method, 'tsreal') ...
                        || strcmpi(proj_method, 'osimagoutput') ...
                        || strcmpi(proj_method, 'osrealoutput')
                    WAV = zeros(size(V, 2), lngthA * size(V, 2));

                    for k = 1:lngthA
                        WAV(:, (k-1)*sizeV+1:k*sizeV) = W' * (sys.A{k} * V);
                    end

                    [Wtmp, ~, ~] = svd(WAV, 'econ');

                    if strcmpi(proj_method, 'osimagoutput') ...
                            || strcmpi(proj_method, 'osrealoutput')
                        Vtmp = Wtmp;
                    end
                end
                
                time_minrel_subspaces = toc(time_minrel_subspaces);
            end
        elseif strcmpi(method, 'strint_linf')
            if strcmpi(proj_method, 'tsimag') ...
                    || strcmpi(proj_method, 'tsreal') ...
                    || strcmpi(proj_method, 'osimaginput') ...
                    || strcmpi(proj_method, 'osrealinput')
                V = tmp.V;
            end
            if strcmpi(proj_method, 'osimaginput') ...
                    || strcmpi(proj_method, 'osrealinput')
                W = V;
            end
            if strcmpi(proj_method, 'tsimag') ...
                    || strcmpi(proj_method, 'tsreal') ...
                    || strcmpi(proj_method, 'osimagoutput') ...
                    || strcmpi(proj_method, 'osrealoutput')
                W = tmp.W;
            end
            if strcmpi(proj_method, 'osimagoutput') ...
                    || strcmpi(proj_method, 'osrealoutput')
                V = W;
            end
            tf = tmp.tf;
        end
    else
        error('Presampling data not found')
    end
end
fprintf(1, '\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ROM.                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, 'Compute ROM.\n');
fprintf(1, '------------\n');
for r=rs
    fprintf(1, 'r = %d ... ', r);
    time_step = tic;
    
    time_mor = tic;
    if strcmp(method, 'strint_avg')
        [Ar, br, cr, result.ctime] = strint_avg(sys, V, W, r, proj_method);
    elseif strcmp(method, 'strint_equi')
        if contains(proj_method, 'imag')
            ns = r;
        elseif contains(proj_method, 'real')
            if mod(r,2)
                fprintf(1, 'skipped.\n')
                continue
            end
            ns = r/2;
        end
        [Ar, br, cr, w, result.ctime] = strint_equi(sys, wmin, wmax, w_pre, ns, proj_method);
    elseif strcmp(method, 'strint_linf')
        [Ar, br, cr, w, result.ctime] = strint_linf(sys, w, V, W, tf, r, proj_method);
    elseif strcmp(method, 'minrel')
        [Ar, br, cr, result.ctime] = minrel(sys, V, W, Vtmp, Wtmp, r, proj_method);
        result.ctime.minrel_dominant_subspaces = time_minrel_subspaces;
    elseif strcmp(method, 'sobt')
        [Ar, br, cr, w, result.ctime] = sobt(sys, r, proj_method, bench);
    end
    result.ctime.mor = toc(time_mor);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE RES.                                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if startsWith(sys.name, 'plate_48')
        n_nodes = full(sum(sum(sys.c)));
    end
    s   = sys.s;

    time_evaluate_rom = tic;
    res = zeros(1, length(s));
    for k=1:length(s)
        Atmp = sys.fA{1}(s(k)) * Ar{1};
        for j=2:length(Ar)
            Atmp = Atmp + sys.fA{j}(s(k)) * Ar{j};
        end

        if isfield(sys, 'fb')
            btmp = sys.fb(s(k)) * br;
        else
            btmp = br;
        end

        if startsWith(sys.name, 'plate_48') && ~contains(sys.name, 'single')
            tmp = cr * (Atmp \ btmp);
            res(k) = sqrt((tmp'*tmp)/n_nodes);
        else
            res(k) = cr * (Atmp \ btmp);
        end
    end
    result.ctime.evaluate_rom = toc(time_evaluate_rom);


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE RESULTS.                                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist('results', 'dir') ~= 7
        mkdir('results')
    end
    
    result.ns = ns;
    result.r = r;
    result.w = w;
    result.res = res;
    
    if isempty(pre_method)
        fname = ['results/' bench '_' method '_' proj_method '.mat'];
    else
        if exist('r0','var') == 1
            fname = ['results/' bench '_' method '_' pre_method '_' ...
                proj_method '_r_' num2str(r0) '.mat'];
        else
            fname = ['results/' bench '_' method '_' pre_method '_' proj_method '.mat'];
        end
    end

    mfile = matfile(fname, 'Writable', true);
    mfile.result_data(1, r) = result;
    
    fprintf(1, 'Completed in %.3f s at %s\n', toc(time_step), datetime('now'));

end

fprintf(1, '\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINISHED SCRIPT.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '\tCompleted in %.3fs at %s\n\n', toc(ctime_all), datetime('now'));
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

end
