function presample_aaa_arnoldi(sys, r, wmin, wmax, ns, nsmin, nsmax, w_pre, side)
%PRESAMPLE_AAA_ARNOLDI Presample interpolation bases for AAAa
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
%   r - size of the Krylov subspace
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

if nargin == 5
    nsmin = 1;
    nsmax = ns;
    w_pre = [];
end
if nargin == 7
    w_pre = [];
end
if nargin ~= 9
    side = 'ts';
end

fprintf('Presampling process PID=%d running on node %s with %d threads.\n',...
    feature('getpid'),...
    strtrim(evalc('system(''hostname'');')),...
    maxNumCompThreads);

%% Sampling points and sizes
w  = linspace(wmin, wmax, ns - length(w_pre)) * 1i;
w  = sort([w 1i*w_pre]);
n = size(sys.A{1},1);

% check folder
if exist('results', 'dir') ~= 7; mkdir('results');  end

if nsmin == 1 && nsmax == ns
    fname = ['presampling/presampling_aaaa_' sys.name '.mat'];
else
    fname = ['presampling/presampling_aaaa_' sys.name '_' num2str(nsmin) ...
        '_' num2str(nsmax) '.mat'];
end

fprintf(1, 'Start presampling %s\n\n', sys.name);

%% Compute bases
tf = zeros(1,1,ns);
time_all_bases = tic;
for k = nsmin:nsmax
    fprintf(1, 'Presampling step %3d / %3d ... ', k, nsmax);
    time_basis = tic;
    
    s = w(k);
    
    % AAA
    nfuns = size(sys.symfuns.dsfuns,2);
    A_aaa = cell(1,nfuns);
    b_aaa = cell(1,nfuns);
    c_aaa = cell(1,nfuns);

    for ii=1:nfuns
        % skip constant functions
        if logical(sys.symfuns.dsfuns{2,ii} == 0); continue; end

        % sampling for AAA
        Z = 1i*linspace(wmin, wmax, 100);
        F = double(sys.symfuns.dsfuns{1,ii}(Z));

        % AAA in matrix form
        [A_aaa{ii},b_aaa{ii},c_aaa{ii}] = getAAAapproximate(F,Z,s);
    end
    
    % matrices for Arnoldi
    At = cell(1,3);
    for ii=1:length(At)
        At{ii} = sparse(n,n);
    end

    for ii=1:nfuns
        At{1} = At{1} + double(sys.symfuns.dsfuns{1,ii}(s))*sys.A{ii};
        if logical(sys.symfuns.dsfuns{2,ii} ~= 0)
            for jj=1:2
                At{1+jj} = At{1+jj} + c_aaa{ii}*A_aaa{ii}^jj*b_aaa{ii}*sys.A{ii};
            end
        end
    end
    
    if isfield(sys.symfuns, 'dsfunb')
        b = double(sys.symfuns.dsfunb{1}(s)) * sys.b;
    else
        b = sys.b;
    end
    
    % directly save V and W to file (tf not possible, because a 1x1xK array
    % is converted to a 1xK array and the indexing throws an error)
    mfile = matfile(fname, 'Writable', true);
    
    % Compute bases and transfer functions
    % Accumulate inputs.
    if strcmpi(side, 'ts') || strcmpi(side, 'input')
        Afuns = cell(1);
        [L, U, P, Q, R] = lu(At{1});
        v0k = Q*(U\(L\(P*(R\b))));
        for ii=2:length(At)
            Afuns{ii-1} = @(x) -Q*(U\(L\(P*(R\(At{ii}*x)))));
        end
        x = TOAR(Afuns{1},Afuns{2},n,r,v0k);
        
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
        Bfuns = cell(1);
        [L, U, P, Q, R] = lu(At{1}.');
        w0k = Q*(U\(L\(P*(R\sys.c'))));
        for ii=2:length(At)
            Bfuns{ii-1} = @(x) -Q*(U\(L\(P*(R\(At{ii}'*x)))));
        end
        y = TOAR(Bfuns{1},Bfuns{2},n,r,w0k);
        
        idx = (size(y,2) * (k - 1) + 1):(size(y,2) * k);
        mfile.W(1:n, idx) = y;
    end
    
    fprintf(1, 'Completed in %.3f s at %s\n', toc(time_basis), datetime('now'));
end
ctime = toc(time_all_bases);

save(fname, 'tf', 'w', 'ctime' ,'-append');

end

