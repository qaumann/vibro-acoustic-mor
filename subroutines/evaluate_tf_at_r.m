% Evaluate the transfer function for a reduced order model with defined size
% r and saves it to csv, fig, and png

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

clear; close all

%% add paths
addpath(genpath('models'))
addpath('utilities')

%% define parameters and output paths
r = 220;

% bench = 'radiation';
% bench = 'transmission';
% bench = 'poroacoustic';
% bench = 'plate_48_rayleigh';
% bench = 'plate_48_hysteretic';
bench = 'plate_48_rayleigh_single';

dir_results = 'results';
dir_graphics= dir_results;
dir_data = dir_results;

save_graphics_flag = true;
save_data_flag = true;

sys = load_model(bench);
if size(sys.res,1) > 1
    sys.res = sys.res(5,:);
end

%% plots
methods = ["strint_linf_aaaa"];
proj_methods = ["tsimag", "tsreal", "osimaginput", "osrealinput", ...
   "osimagoutput", "osrealoutput"];

for m = methods
    for pm = proj_methods
        load([dir_results '/' bench '_' m{:} '_' pm{:} '.mat'], 'result_data');
        
        fig = figure(1);
        set(fig,'defaulttextinterpreter','latex')
        subplot(2,1,1)
        semilogy(abs(sys.s)/2/pi, abs(sys.res), ...
            abs(sys.s)/2/pi, abs(result_data(end).res(1:length(sys.res))));
        xlabel('Frequency $\omega$ (rad/sec)');
        ylabel('Magnitude');
        legend({'Full model', ['$r=' num2str(result_data(r).r) '$']}, 'interpreter', 'latex');
        title(['Transfer functions -- ' replace(m{:}, '_', ' ')]);

        subplot(2, 1, 2);
        if xor(all(imag(sys.res) == 0), all(imag(result_data(r).res) == 0))
            hinferr = abs(abs(sys.res) - ...
                abs(result_data(r).res(1:length(sys.res)))) ...
                ./ abs(sys.res);
            semilogy(abs(sys.s)/2/pi, hinferr);
        else
            hinferr = abs(sys.res - result_data(r).res(1:length(sys.res))) ...
                ./ abs(sys.res);
            semilogy(abs(sys.s)/2/pi, hinferr);
        end

        xlabel('Frequency $\omega$ (rad/sec)');
        ylabel('Magnitude');
        title('Relative errors');

        if save_graphics_flag == true
            saveas(gcf, [dir_graphics '/' bench '_tf_' m{:} '_' pm{:} '_r_' num2str(r) '.png'], 'png');
            saveas(gcf, [dir_graphics '/' bench '_tf_' m{:} '_' pm{:} '_r_' num2str(r) '.fig'], 'fig');
        end
        close(1)

        if save_data_flag == true
            data = [abs(sys.s)/2/pi; abs(sys.res); ...
                abs(result_data(end).res(1:length(sys.res))); hinferr];
            qcsv([dir_data '/' bench '_tf_' m{:} '_' pm{:} '_r_' num2str(r) '.csv'], transpose(data), ...
                'header', 's,res,res_r,hinferr' );
        end
    end
end
