function [Ar, br, cr, w, ctime] = strint_linf(sys, w, V, W, tf, r, method)
%STRINT_LINF Compute structured interpolation with L-infinity points.
%
% SYNTAX:
%   [Ar, br, cr, w, ctime] = STRINT_LINF(sys, w, V, W, tf, r, method)
%
% DESCRIPTION: Requires presampled data (see presample.m)
%
% INPUT:
%   sys - system data (see load_model.m)
%   w - expansion frequencies
%   V, W - presampling bases (see presample.m)
%   tf - precomputed transfer function data (see presample.m)
%   r - size of the reduced model
%   method - type of projection space for interpolation
%              'tsimag'       - two-sided projection allowing imag. parts
%              'tsreal'       - two-sided projection only real parts
%              'osimaginput'  - one-sided projection imag. input parts
%              'osrealinput'  - one-sided projection real input parts
%              'osimagoutput' - one-sided projection imag. output parts
%              'osrealoutput' - one-sided projection real output parts
%
% OUTPUT:
%   Ar, br, cr - reduced order model
%   w - expansion frequencies
%   ctime - struct with computation times

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Sizes
n  = size(sys.A{1}, 1);
nA = length(sys.A);
ns = length(w);
rs = size(V,2) / length(w);

if isfield(sys, 'fb')
    nb = length(sys.fb);
else
    nb = 0;
end

if isfield(sys, 'fc')
    nc = length(sys.fc);
else
    nc = 0;
end
if startsWith(sys.name, 'plate_48')
    n_nodes = full(sum(sum(sys.c)));
end

ctime = struct();

% Compute final bases.
k            = 1;
idx          = zeros(1, r);
Vtmp         = zeros(n, 0);
Wtmp         = zeros(n, 0);
[Ar, br, cr] = compute_rom(sys, nb, nc, Vtmp, Wtmp);
while k <= ceil(r/rs)
    time_errsys = tic;
    errsys = zeros(1, ns);
    for j = 1:ns
        if ismember(j, idx)
            errsys(j) = 0;
        else
            s = w(j);
            
            % Accumulate system matrices.
            Atmp = sys.fA{1}(s) * Ar{1};
            for i = 2:nA
                Atmp = Atmp + sys.fA{i}(s) * Ar{i};
            end
            
            if nb == 0
                btmp = br;
            elseif nb == 1
                btmp = sys.fb(s) * br;
            else
                btmp = sys.fb{1}(s) * br{1};
                for i = 2:nb
                    btmp = btmp + sys.fb{i}(s) * br{i};
                end
            end
            
            if nc == 0
                ctmp = cr;
            elseif nc == 1
                ctmp = sys.fc(s) * cr;
            else
                ctmp = sys.fc{1}(s) * cr{1};
                for i = 2:nc
                    ctmp = ctmp + sys.fc{i}(s) * cr{i};
                end
            end
            
            % Compute error system.
            if startsWith(sys.name, 'plate_48') && ~contains(sys.name, 'single')
                tmp = ctmp * (Atmp \ btmp);
                errsys(j) = max(svd(tf(:, :, j) - sqrt((tmp'*tmp)/n_nodes)));
            else
                errsys(j) = max(svd(tf(:, :, j) - ctmp * (Atmp \ btmp)));
            end
        end
    end
    ctime.errsys(k) = toc(time_errsys);
    
    [~, tmp] = max(errsys);
    idx(k)   = tmp;
    
    select   = idx(1:k);
    select   = select(not(select == 0));
    new_select = [];
    for ii=select
        new_select = [new_select ii*rs-(rs-1):ii*rs]; %#ok
    end
    select = new_select;
    
    % Compute next ROM.
    time_qr = tic;
    if strcmpi(method, 'tsimag') || strcmpi(method, 'osimaginput')
        [Vtmp, ~, ~] = qr(V(:, select), 0);
    elseif strcmpi(method, 'tsreal') || strcmpi(method, 'osrealinput')
        [Vtmp, ~, ~] = qr([real(V(:, select)), imag(V(:, select))], 0);
    end
    if strcmpi(method, 'osimaginput') || strcmpi(method, 'osrealinput')
        Wtmp = Vtmp;
    end
    
    if strcmpi(method, 'tsimag') || strcmpi(method, 'osimagoutput')
        [Wtmp, ~, ~] = qr(W(:, select), 0);
    elseif strcmpi(method, 'tsreal') || strcmpi(method, 'osrealoutput')
        [Wtmp, ~, ~] = qr([real(W(:, select)), imag(W(:, select))], 0);
    end
    if strcmpi(method, 'osimagoutput') || strcmpi(method, 'osrealoutput')
        Vtmp = Wtmp;
    end
    ctime.qr(k) = toc(time_qr);
    
    time_projection = tic;
    [Ar, br, cr] = compute_rom(sys, nb, nc, Vtmp, Wtmp);
    ctime.projection(k) = toc(time_projection);
    
    if strcmpi(method, 'tsreal') ...
            || strcmpi(method, 'osrealinput') ...
            || strcmpi(method, 'tsreal') ...
            || strcmpi(method, 'osrealoutput')
        k = k + 2;
    else
        k = k + 1;
    end
end

end

% Compute reduced-order model.
function [Ar, br, cr] = compute_rom(sys, nb, nc, V, W)

Ar = cellfun(@(c) W' * (c * V), sys.A, 'UniformOutput', 0);

if nb > 1
    br = cellfun(@(c) W' * c, sys.b, 'UniformOutput', 0);
else
    br = W' * sys.b;
end

if nc > 1
    cr = cellfun(@(c) c * V, sys.c, 'UniformOutput', 0);
else
    cr = sys.c * V;
end

end
