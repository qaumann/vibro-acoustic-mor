function [Abar,bbar,c] = getAAAapproximate(F,Z,s0)
%GETAAAAPROXIMATE compute AAA interpolation in shifted matrix form

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

[w,l,z,~] = aaa(F,Z,1e-14,50);
[c,A,E,b] = aaa2mat(w,l,z);
Abar = -(A+s0*E)\E;
bbar = (A+s0*E)\b;
end
