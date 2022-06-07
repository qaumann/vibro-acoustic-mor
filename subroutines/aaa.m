function [w,f,z,errvec] = aaa(F,Z,tol,mmax)
%AAA compute AAA interpolation
%
%   Reference: Nakatsukasa, Sète, Trefethen. (2018). The AAA algorithm for 
%       rational approximation. SIAM Journal on Scientific Computing,
%       40(3), A1494-A1522.

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

M = length(Z);
if nargin<3, tol=1e-13; end
if nargin<4, mmax=100; end

%column vectors
Z=Z(:);
F=F(:);

SF = spdiags(F,0,M,M);
J = 1:M;
R = mean(F);
z=[];f=[];C=[];errvec=[];

for m=1:mmax
    %select next support point
    [~,j] = max(abs(F-R));
    z = [z; Z(j)];  %#ok
    f = [f; F(j)];  %#ok
    J(J==j) = [];
    
    %Cauchy matrix
    C = [C 1./(Z-Z(j))]; %#ok
    
    %right scaling matrix
    Sf = diag(f);
    
    %Loewner matrix
    A = SF*C - C*Sf;
    [~,~,V] = svd(A(J,:),0);
    w = V(:,m);
    N = C*(w.*f);
    D = C*w;
    R = F;
    R(J) = N(J)./D(J);
    err = norm(F-R,inf);
    errvec = [errvec; err]; %#ok
    if err <= tol*norm(F,inf); break; end
end

end

