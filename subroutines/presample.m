function presample(sys, wmin, wmax, ns, nsmin, nsmax, w_pre, side)
%PRESAMPLE Presample interpolation bases for MINREL, STRING_AVG, STRINT_LINF
%
% SYNTAX: [V, W, tf] = presample(sys, wmin, wmax, ns, nsmin, nsmax)
%
% DESCRIPTION: Computes complex-valued left and right interpolation bases
%   and the transfer function for all interpolation points. Real bases have
%   to be computed from V, W in separate algorithms.
%   The dominant subspaces Vtmp, Wtmp for minrel have to be computed in a
%   separate function.
%
% INPUT:
%   sys - full order system
%   wmin, wmax - bounds for interpolation points
%   ns - number of sampling points
%
% OPTIONAL:
%   nsmin, nsmax - indices where to start and to end if only slices of the
%       basis are computed
%   w_pre - predefined sampling points
%   side - which side to be considered (default: twosided). Can be 'input' 
%       or 'output'
%
% OUTPUT:
%   V, W - left and right bases
%   tf - transfer function evaluations

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin == 4
    nsmin = 1;
    nsmax = ns;
    w_pre = [];
end
if nargin == 6
    w_pre = [];
end
if nargin ~= 8
    side = 'ts';
end

fprintf('Presampling process PID=%d running on node %s with %d threads.\n',...
    feature('getpid'),...
    strtrim(evalc('system(''hostname'');')),...
    maxNumCompThreads);

%% Sampling points and sizes
w  = linspace(wmin, wmax, ns - length(w_pre)) * 1i;
w  = sort([w 1i*w_pre]);
n  = size(sys.A{1}, 1);
nA = length(sys.A);

if isfield(sys, 'fb')
    nb = length(sys.fb);
else
    nb = 0;
end

if isa(sys.b, 'cell')
    m = size(sys.b{1}, 2);
else
    m = size(sys.b, 2);
end

if isfield(sys, 'fc')
    nc = length(sys.fc);
else
    nc = 0;
end

if isa(sys.c, 'cell')
    p = size(sys.c{1}, 1);
else
    p = size(sys.c, 1);
end

tf = zeros(m, p, ns);

% check folder
if exist('presampling', 'dir') ~= 7; mkdir('presampling');  end

if nsmin == 1 && nsmax == ns
    fname = ['presampling/presampling_' sys.name '.mat'];
else
    fname = ['presampling/presampling_' sys.name '_' num2str(nsmin) ...
        '_' num2str(nsmax) '.mat'];
end

fprintf(1, 'Start presampling %s\n\n', sys.name);

%% Compute bases
% time_assembly = zeros(
time_all_bases = tic;
for k = nsmin:nsmax
    fprintf(1, 'Presampling step %3d / %3d ... ', k, nsmax);
    time_basis = tic;
    time_assembly = tic;
    
    s = w(k);
    
    % Accumulate system matrices.
    Atmp = sys.fA{1}(s) * sys.A{1};
    for j = 2:nA
        Atmp = Atmp + sys.fA{j}(s) * sys.A{j};
    end

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
    
    time_assembly = toc(time_assembly);
    
    % directly save V and W to file (tf not possible, because a 1x1xK array
    % is converted to a 1xK array and the indexing throws an error)
    mfile = matfile(fname, 'Writable', true);
    
    % Compute bases and transfer functions
    time_solve = tic;
    if strcmpi(side, 'ts') || strcmpi(side, 'input')
        x = full(Atmp \ btmp);
        mfile.V(1:n, k) = x;
        tf(:, :, k)     = ctmp * x;
    end
    if strcmpi(side, 'ts') || strcmpi(side, 'output')
        y = full(Atmp' \ ctmp');
        mfile.W(1:n, k) = y;
        tf(:, :, k)     = y' * btmp;
    end
    mfile.ctime_solve(1,k) = toc(time_solve);
    mfile.ctime_assembly(1,k) = time_assembly;
    mfile.ctime_basis(1,k) = toc(time_basis);
    
    fprintf(1, 'Completed in %.3f s at %s\n', mfile.ctime_basis(1,k), datetime('now'));
end
ctime = toc(time_all_bases);

save(fname, 'tf', 'w', 'ctime' ,'-append');

end

