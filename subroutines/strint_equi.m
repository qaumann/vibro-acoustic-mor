function [Ar, br, cr, w, ctime] = strint_equi(sys, wmin, wmax, w_pre, ns, method)
%STRINT_EQUI Compute structured interpolation with equidistant points.
%
% SYNTAX:
%   [Ar, br, cr] = STRINT_EQUI(sys, wmin, wmax, ns, method)
%
% DESCRIPTION:
%
% INPUT:
%   sys - system data (see load_model.m)
%   wmin, wmax - min and max expansion frequency
%   w_pre - predefined expansion frequencies (can be empty)
%   ns - number of equidistant expansion points
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

if length(w_pre) > ns
    error('Cannot prescribe more than ns points')
end

% Interpolation points.
if ns==1
    w = 1i * (wmin + wmax)/2;
else
    w = linspace(wmin, wmax, ns - length(w_pre)) * 1i;
    w = sort([w 1i*w_pre]);
end
n  = size(sys.A{1}, 1);
nA = length(sys.A);

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

ctime = struct();

% Interpolation bases sizes.
switch lower(method)
    case 'tsimag'
        V = zeros(n, ns);
        W = zeros(n, ns);
        
    case 'tsreal'
        V = zeros(n, 2 * ns);
        W = zeros(n, 2 * ns);
        
    case 'osimaginput'
        V = zeros(n, ns);
        
    case 'osrealinput'
        V = zeros(n, 2 * ns);
        
    case 'osimagoutput'
        W = zeros(n, ns);
        
    case 'osrealoutput'
        W = zeros(n, 2 * ns);
        
    otherwise
        error('Requested projection method not implemented!');
end

% Compute bases
for k = 1:ns
    fprintf(1, 'Step %3d / %3d\n', k, ns);
    time_basis = tic;
    time_assemble = 0;
    timer_assemble = tic;
    time_solve = 0;
    
    s = w(k);
    
    % Accumulate system matrix.
    Atmp = sys.fA{1}(s) * sys.A{1};
    for j = 2:nA
        Atmp = Atmp + sys.fA{j}(s) * sys.A{j};
    end
    
    time_assemble = time_assemble + toc(timer_assemble);
    
    % Accumulate inputs.
    if strcmpi(method, 'tsimag') ...
            || strcmpi(method, 'tsreal') ...
            || strcmpi(method, 'osimaginput') ...
            || strcmpi(method, 'osrealinput')
        timer_assemble = tic;
        if nb == 0
            btmp = sys.b;
        elseif nb == 1
            btmp = sys.fb(s) * sys.b;
        else
            btmp = sys.fb{1}(s) * sys.b{1};
            for j = 2:nb
                btmp = btmp + sys.fb{j}(s) * sys.b{j};
            end
        end
        time_assemble = time_assemble + toc(timer_assemble);
        
        timer_solve = tic;
        x = full(Atmp \ btmp);
        time_solve = time_solve + toc(timer_solve);
        
        if strcmpi(method, 'tsimag') || strcmpi(method, 'osimaginput')
            V(:, k) = x;
        else
            V(:, 2*k-1:2*k) = [real(x), imag(x)];
        end
    end
    
    % Accumulate outputs.
    if strcmpi(method, 'tsimag') ...
            || strcmpi(method, 'tsreal') ...
            || strcmpi(method, 'osimagoutput') ...
            || strcmpi(method, 'osrealoutput')
        timer_assemble = tic;
        if nc == 0
            ctmp = sys.c;
        elseif nc == 1
            ctmp = sys.fc(s) * sys.c;
        else
            ctmp = sys.fc{1}(s) * sys.c{1};
            for j = 2:nc
                ctmp = ctmp + sys.fc{j}(s) * sys.c{j};
            end
        end
        time_assemble = time_assemble + toc(timer_assemble);
        
        timer_solve = tic;
        y = full(Atmp' \ ctmp');
        time_solve = time_solve + toc(timer_solve);
        
        if strcmpi(method, 'tsimag' ) || strcmpi(method, 'osimagoutput')
            W(:, k) = y;
        else
            W(:, 2*k-1:2*k) = [real(y), imag(y)];
        end
    end
    
    ctime.time_assemble(k) = time_assemble;
    ctime.time_solve(k) = time_solve;
    ctime.time_basis(k) = toc(time_basis);
end

time_qr = tic;
if strcmpi(method, 'tsimag') ...
        || strcmpi(method, 'tsreal') ...
        || strcmpi(method, 'osimaginput') ...
        || strcmpi(method, 'osrealinput')
    [V, ~, ~] = qr(V, 0);
    
    if strcmpi(method, 'osimaginput') || strcmpi(method, 'osrealinput')
        W = V;
    end
end

if strcmpi(method, 'tsimag') ...
        || strcmpi(method, 'tsreal') ...
        || strcmpi(method, 'osimagoutput') ...
        || strcmpi(method, 'osrealoutput')
    [W, ~, ~] = qr(W, 0);
    
    if strcmpi(method, 'osimagoutput') || strcmpi(method, 'osrealoutput')
        V = W;
    end
end
ctime.qr = toc(time_qr);

% Compute reduced-order model.
time_projection = tic;
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
ctime.projection = toc(time_projection);

end
