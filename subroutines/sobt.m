function [Ar, br, cr, w] = sobt(sys, r, method, bench)
%SOBT Compute second-order balanced truncation.
%
% method: balancing formulas 1--8, osinput, osoutput

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

addpath(genpath('software'))

assert(length(sys.fA) <= 3, ...
    'SOBT cannot be used for systems different from second-order form!');

n = size(sys.A{1}, 1);

% add zero A{2} (damping) matrix if required
if length(sys.A) == 2
    sys.A{3} = sys.A{2};
    sys.A{2} = sparse(n,n);
    sys.fA{3} = sys.fA{2};
    sys.fA{2} = sys.fA{1};
end

% Compute controllability Gramian.
fprintf(1, 'Compute controllability Gramian.\n');
if not(strcmpi(method, 'osoutput'))
    if exist(['presampling/' bench '_gram_c.mat'], 'file')
        tmp = load(['presampling/' bench '_gram_c.mat']);
        ZC  = tmp.ZC;
    else
        n           = size(sys.A{1}, 1);
        opts        = struct('norm', 'fro');
        opts.adi    = struct( ...
            'maxiter'     , 1500, ...
            'res_tol'     , 1e-10, ...
            'rel_diff_tol', 0, ...
            'computeZ'    , 1, ...
            'info'        , 1);
        opts.shifts = struct( ...
            'info'       , 0, ...
            'num_desired', 10, ...
            'method'     , 'projection');
        oper        = operatormanager('default');
        
        eqn  = struct( ...
            'A_', [sparse(n, n), speye(n); -sys.A{1}, -sys.A{2}], ...
            'B' , [zeros(size(sys.b)); sys.b], ...
            'C' , [sys.c, zeros(size(sys.c))], ...
            'E_', [speye(n), sparse(n, n); sparse(n, n), sys.A{3}]);
        eqn.type  = 'N';
        eqn.haveE = 1;
        
        opts.shifts.p = mess_para(eqn, opts, oper);
        out           = mess_lradi(eqn, opts, oper);
        ZC            = out.Z;
        
        if exist('presampling', 'dir') ~= 7; mkdir('presampling');  end
        save(['presampling/' bench '_gram_c.mat'], 'ZC', '-v7.3');
    end
end

% Compute observability Gramian.
fprintf(1, 'Compute observability Gramian.\n');
if not(strcmpi(method, 'osinput'))
    if exist(['presampling/' bench '_gram_o.mat'], 'file')
        tmp = load(['presampling/' bench '_gram_o.mat']);
        ZO  = tmp.ZO;
    else
        n           = size(sys.A{1}, 1);
        opts        = struct('norm', 'fro');
        opts.adi    = struct( ...
            'maxiter'     , 1500, ...
            'res_tol'     , 1e-10, ...
            'rel_diff_tol', 0, ...
            'computeZ'    , 1, ...
            'info'        , 1);
        opts.shifts = struct( ...
            'info'       , 0, ...
            'num_desired', 10, ...
            'method'     , 'projection');
        oper        = operatormanager('default');
        
        eqn  = struct( ...
            'A_', [sparse(n, n), speye(n); -sys.A{1}, -sys.A{2}], ...
            'B' , [zeros(size(sys.b)); sys.b], ...
            'C' , [sys.c, zeros(size(sys.c))], ...
            'E_', [speye(n), sparse(n, n); sparse(n, n), sys.A{3}]);
        eqn.type  = 'T';
        eqn.haveE = 1;
        
        opts.shifts.p = mess_para(eqn, opts, oper);
        out           = mess_lradi(eqn, opts, oper);
        ZO            = out.Z;
        
        if exist('presampling', 'dir') ~= 7; mkdir('presampling');  end
        save(['presampling/' bench '_gram_o.mat'], 'ZO', '-v7.3');
    end
end

% Balancing or dominant subspaces.
fprintf(1, 'Compute ROM.\n');
w = [];
if strcmpi(method, 'osinput')
    [V, ~, ~] = qr(ZC, 0);
    [V, ~, ~] = qr(V(1:n, 1:r), 0);
    W         = V;
    
    Ar = cellfun(@(c) W' * (c * V), sys.A, 'UniformOutput', 0);
    br = W' * sys.b;
    cr = sys.c * V;
elseif strcmpi(method, 'osoutput')
    [W, ~, ~] = qr(ZO, 0);
    [W, ~, ~] = qr(W(1:n, 1:r), 0);
    V         = W;
    
    Ar = cellfun(@(c) W' * (c * V), sys.A, 'UniformOutput', 0);
    br = W' * sys.b;
    cr = sys.c * V;
else
    tmp  = struct('M', sys.A{3}, 'E', sys.A{2}, 'K', sys.A{1}, ...
        'Bu', sys.b, 'Cp', sys.c);
    opts = struct( ...
        'BalanceType'     , method, ...
        'Method'          , 'sr', ...
        'OrderComputation', 'order', ...
        'Order'           , r);
    tmp  = btred_so(tmp, speye(n), ZC, ZO, opts);
    
    Ar = {tmp.K, tmp.E, tmp.M};
    br = tmp.Bu;
    cr = tmp.Cp;
end
