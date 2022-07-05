clear; close all

experiments = {'plate_48_rayleigh_single_strint_equi_tsimag'};
experiments = {'poroacoustic_strint_equi_osimaginput', ...
    'poroacoustic_strint_equi_osimagoutput', ...
    'poroacoustic_strint_equi_tsimag', ...
    'poroacoustic_strint_equi_tsreal', ...
    'radiation_strint_equi_osimaginput', ...
    'radiation_strint_equi_osimagoutput', ...
    'radiation_strint_equi_tsimag', ...
    'transmission_strint_equi_tsimag'};

if exist(['results' filesep 'backup'],'dir')  ~= 7
    mkdir(['results' filesep 'backup'])
end
files = dir('results');

for exp = experiments
    copyfile(['results' filesep exp{:} '.mat'], ['results' filesep 'backup'])
    
    load(['results' filesep exp{:} '.mat'], 'result_data');
    inds = startsWith({files.name}, [exp{:} '_']);
    
    partial_files = files(inds);
    
    for ii = 1:length(partial_files)
        tmp = load(['results' filesep partial_files(ii).name]);
        ninds = [tmp.result_data(:).r];
        result_data(ninds) = tmp.result_data(ninds);
        movefile(['results' filesep partial_files(ii).name], ['results' filesep 'backup'])
    end
    
    save(['results' filesep exp{:} '.mat'], 'result_data');
end
