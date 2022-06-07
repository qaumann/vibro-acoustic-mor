function presample_structure_preserving(sys, wmin, wmax, ns, nsmin, nsmax, w_pre, side)
%PRESAMPLE_STRUCTURE_PRESERVING Presample interpolation bases for structure preserving
%
% SYNTAX: presample(sys, wmin, wmax, ns, nsmin, nsmax)
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
% OUTPUT: The output is directly saved to a MAT file
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

addpath('subroutines')

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
n = size(sys.A{1}, 1);
nfuns = size(sys.symfuns.dsfuns,2);

r = size(sys.symfuns.dsfuns,1);
K = cell(1,r);

% check folder
if exist('presampling', 'dir') ~= 7; mkdir('presampling');  end

if nsmin == 1 && nsmax == ns
    fname = ['presampling/presampling_strprs_' sys.name '.mat'];
else
    fname = ['presampling/presampling_strprs_' sys.name '_' num2str(nsmin) ...
        '_' num2str(nsmax) '.mat'];
end

fprintf(1, 'Start presampling %s\n\n', sys.name);

%% Compute bases
tf = zeros(1,1,ns);
time_all_bases = tic;
for k = nsmin:nsmax
    fprintf(1, 'Presampling step %3d / %3d ... ', k, nsmax);
    time_basis = tic;
    time_assemble = tic;
    
    s = w(k);
    
    % matrices for recursion
    for ii=1:length(K)
        K{ii} = sparse(n,n);
    end
    
    for ii=1:nfuns
        K{1} = K{1} + double(sys.symfuns.dsfuns{1,ii}(s))*sys.A{ii};
        for jj = 1 : r-1
            factor = double(sys.symfuns.dsfuns{jj+1,ii}(s)) / factorial(jj);
            if abs(factor) > eps
                K{1+jj} = K{1+jj} + factor * sys.A{ii};
            end
        end
        
    end
    
    time_assemble = toc(time_assemble);
    
    time_decompose = tic;
    K{1} = decomposition(K{1});
    time_decompose = toc(time_decompose);
    
    % directly save V and W to file (tf not possible, because a 1x1xK array
    % is converted to a 1xK array and the indexing throws an error)
    mfile = matfile(fname, 'Writable', true);
    
    % Compute bases and transfer functions
    time_solve = tic;
    
    % Accumulate inputs.
    if strcmpi(side, 'ts') || strcmpi(side, 'input')
        x = zeros(n,r);
        if isfield(sys.symfuns, 'dsfunb')
            x(:,1) = K{1} \ (double(sys.symfuns.dsfunb{1}(s)) .* sys.b);
        else
            x(:,1) = K{1} \ sys.b;
        end
        for jj=2:r
            if isfield(sys.symfuns, 'dsfunb')
                factorb = double(sys.symfuns.dsfunb{jj}(s)) / factorial(jj-1);
                tmpb = factorb * sys.b;
            else
                tmpb = zeros(n,1);
            end
            for kk=1:jj-1
                tmpb = tmpb - K{kk+1} * x(:,jj-kk);
            end
            x(:,jj) = K{1} \ tmpb;
        end
        idx = (size(x,2) * (k - 1) + 1):(size(x,2) * k);
        mfile.V(1:n, idx) = x;
    
        % compute transfer function
        if startsWith(sys.name, 'plate') && ~contains(sys.name, 'single')
            n_nodes = full(sum(sum(sys.c)));
            tmp = sys.c * x(:,1);
            tf(:, :, k) = sqrt((tmp'*tmp)/n_nodes);
        else
            tf(:, :, k) = sys.c * x(:,1);
        end
    end
    
    % Accumulate outputs.
    if strcmpi(side, 'ts') || strcmpi(side, 'output')
        y = zeros(n,r);
        y(:,1) = K{1}' \ transpose(sys.c);
        for jj=2:r
            tmpc = zeros(n,1);
            for kk=1:jj-1
                tmpc = tmpc - K{kk+1}' * y(:,jj-kk);
            end
            y(:,jj) = K{1}' \ tmpc;
        end
        
        idx = (size(y,2) * (k - 1) + 1):(size(y,2) * k);
        mfile.W(1:n, idx) = y;
    end
    
    mfile.ctime_solve(1,k) = toc(time_solve);
    mfile.ctime_assemble(1,k) = time_assemble;
    mfile.ctime_aaa(1,k) = time_decompose;
    mfile.ctime_basis(1,k) = toc(time_basis);
    
    fprintf(1, 'Completed in %.3f s at %s\n', mfile.ctime_basis(1,k), datetime('now'));
end
ctime = toc(time_all_bases);

save(fname, 'tf', 'w', 'ctime' ,'-append');

end

