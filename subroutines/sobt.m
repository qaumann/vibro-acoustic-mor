function [Ar, br, cr, w, ctime] = sobt(sys, r, method, bench)
%SOBT Compute second-order balanced truncation.
%
% SYNTAX:
%   [Ar, br, cr, ctime] = STRINT_EQUI(sys, r, method, bench)
%
% DESCRIPTION:
%
% INPUT:
%   sys - system data (see load_model.m)
%   r - size of the reduced model
%   method - type of projection space for interpolation
%              'p'        - position balancing
%              'pm'       - position balancing (diagonalized M)
%              'pv'       - position-velocity balancing
%              'vp'       - velocity-position balancing
%              'vpm'      - velocity-position balancing (diag. M)
%              'v'        - velocity balancing
%              'fv'       - free velocity balancing
%              'so'       - second-order balancing
%              'osinput'  - one-sided projection imag. input parts
%              'osoutput' - one-sided projection imag. output parts
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

ctime = struct();

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
        
        ctime_gram_c = tic;
        opts.shifts.p = mess_para(eqn, opts, oper);
        out           = mess_lradi(eqn, opts, oper);
        ZC            = out.Z;
        ctime_gram_c = toc(ctime_gram_c);
        
        if exist('presampling', 'dir') ~= 7; mkdir('presampling');  end
        save(['presampling/' bench '_gram_c.mat'], 'ZC', 'ctime_gram_c', '-v7.3');
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
        
        ctime_gram_o = tic;
        opts.shifts.p = mess_para(eqn, opts, oper);
        out           = mess_lradi(eqn, opts, oper);
        ZO            = out.Z;
        ctime_gram_o = toc(ctime_gram_o);
        
        if exist('presampling', 'dir') ~= 7; mkdir('presampling');  end
        save(['presampling/' bench '_gram_o.mat'], 'ZO', 'ctime_gram_o', '-v7.3');
    end
end

% Balancing or dominant subspaces.
fprintf(1, 'Compute ROM.\n');
w = [];
if strcmpi(method, 'osinput')
    time_qr = tic;
    [V, ~, ~] = qr(ZC, 0);
    [V, ~, ~] = qr(V(1:n, 1:r), 0);
    W         = V;
    ctime.qr = toc(time_qr);
    
    time_projection = tic;
    Ar = cellfun(@(c) W' * (c * V), sys.A, 'UniformOutput', 0);
    br = W' * sys.b;
    cr = sys.c * V;
    ctime.projection = toc(time_projection);
elseif strcmpi(method, 'osoutput')
    time_qr = tic;
    [W, ~, ~] = qr(ZO, 0);
    [W, ~, ~] = qr(W(1:n, 1:r), 0);
    V         = W;
    ctime.qr = toc(time_qr);
    
    time_projection = tic;
    Ar = cellfun(@(c) W' * (c * V), sys.A, 'UniformOutput', 0);
    br = W' * sys.b;
    cr = sys.c * V;
    ctime.projection = toc(time_projection);
else
    tmp  = struct('M', sys.A{3}, 'E', sys.A{2}, 'K', sys.A{1}, ...
        'Bu', sys.b, 'Cp', sys.c);
    opts = struct( ...
        'BalanceType'     , method, ...
        'Method'          , 'sr', ...
        'OrderComputation', 'order', ...
        'Order'           , r);
    time_bt = tic;
    tmp  = btred_so(tmp, speye(n), ZC, ZO, opts);
    ctime.bt = toc(time_bt);
    
    Ar = {tmp.K, tmp.E, tmp.M};
    br = tmp.Bu;
    cr = tmp.Cp;
end
