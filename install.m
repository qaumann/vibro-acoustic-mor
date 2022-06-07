%% INSTALL.
% Install all software and subroutines.

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

fprintf(1, 'Install software and subroutines.\n');
fprintf(1, '=================================\n');

addpath(genpath('models'));
addpath(genpath('presampling'));
addpath(genpath('results'));
addpath(genpath('slurmscripts'));
addpath(genpath('software'));
addpath(genpath('subroutines'));

fprintf(1, '\n');


%% Finished install.

fprintf(1, 'Finished install.\n');
fprintf(1, '=================\n');
