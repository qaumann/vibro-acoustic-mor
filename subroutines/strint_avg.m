function [Ar, br, cr] = strint_avg(sys, V, W, r, method)
%STRINT_AVG Compute structured interpolation with approximate subspaces.
%
% SYNTAX:
%   [Ar, br, cr] = strint_avg(sys, V, W, r, method)
%
% DESCRIPTION: Requires presampled data (see presample.m)
%
% INPUT:
%   sys - system data (see load_model.m)
%   V, W - presampling bases (see presample.m)
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

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Sizes
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

% Compute final reduction basis
if strcmpi(method, 'tsimag') ...
        || strcmpi(method, 'tsreal') ...
        || strcmpi(method, 'osimaginput') ...
        || strcmpi(method, 'osrealinput')
    [V, ~, ~] = qr(V, 0);
    V         = V(:, 1:r);
    
    if strcmpi(method, 'osimaginput') || strcmpi(method, 'osrealinput')
        W = V;
    end
end

if strcmpi(method, 'tsimag') ...
        || strcmpi(method, 'tsreal') ...
        || strcmpi(method, 'osimagoutput') ...
        || strcmpi(method, 'osrealoutput')
    [W, ~, ~] = qr(W, 0);
    W         = W(:, 1:r);
    
    if strcmpi(method, 'osimagoutput') || strcmpi(method, 'osrealoutput')
        V = W;
    end
end

% Compute reduced-order model.
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
