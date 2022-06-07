function [a,E,F,b] = aaa2mat(w,f,z)
%AAA2MAT convert AAA interpolant to matrix form

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

n = length(w);

a = (f.*w).';
b = eye(n,1);
F = diag([0 -1*ones(1,n-1)]) + diag(ones(1,n-1),-1);
E = diag([0 z(2:end).']) + diag(-z(1:end-1).',-1);
E(1,:) = w.';
end

