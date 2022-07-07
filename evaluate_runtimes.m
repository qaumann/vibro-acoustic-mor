
clear;
close all;

%% define parameters and output paths
benchs = {'plate_48_rayleigh_single', 'plate_48_rayleigh','radiation', 'transmission', 'poroacoustic', ...
     'plate_48_hysteretic'};
proj_methods = ["tsimag", "tsreal", "osimagoutput", "osrealoutput", ...
    "osimaginput", "osrealinput"];

dir_results = 'results';
dir_graphics= 'results_graphics';
dir_data = 'results_data';

save_graphics_flag = true;
save_data_flag = true;

for bb = 1:length(benchs)
    files = dir(dir_results);
    results = struct();
    ind = 1;
    bench = benchs{bb};
    
    for ii=1:length(files)
        if ~startsWith(files(ii).name, bench); continue; end
        if ~endsWith(files(ii).name, '.mat'); continue; end
        %distinguish plate_48_rayleigh from plate_48_rayleigh_single
        if strcmp(bench, 'plate_48_rayleigh') && contains(files(ii).name, 'single')
            continue
        end
        
        name = erase(files(ii).name, [bench + "_", ".mat"]);
        if startsWith(name, 'gram'); continue; end
        if startsWith(name, 'presampling'); continue; end
        if ~(contains(name,'equi') || contains(name,'aaaa') || ...
                contains(name,'strprs') || contains(name,'sobt'))
            tmp = strsplit(name,'_');
            name = [sprintf('%s_',tmp{1:end-1}),'std_',tmp{end}];
        end
        if contains(name,'osinput')
            name = strrep(name,'osinput','osimaginput');
        end
        results(ind).name = name;
        load([files(ii).folder '/' files(ii).name]);
        results(ind).ctime_mor = result_data(end).ctime.mor;
        
        ind = ind + 1;
    end
    
    % add dummy for single-sided bt
    if any(strcmp({results.name},'sobt_osimaginput'))
        results(end+1).name = 'sobt_osrealinput';  %#ok<SAGROW>
        results(end).ctime_mor = 0;
    end
    
    
    [~,idx]=sort({results.name});
    results = results(idx);
    
    ctime_mor = {};
    idx = 1;
    for m = proj_methods
        ind = contains({results.name}, m{:});
        if ~any(ind); continue; end
        this_results = results(ind);
        
        ctime_mor{idx} = [this_results.ctime_mor]'; %#ok
        idx = idx+1;
        
        if save_data_flag == true
                    % typeset names
            this_names = typeset_names({this_results.name}', m{:});
            qcsv([dir_data '/' bench '_ctime_' m{:} '.csv'], ...
                {this_names {this_results.ctime_mor}'}, ...
                'header', 'method,ctime_mor');
        end
    end
    
    bar_names = {results.name};
    for m = proj_methods
        bar_names = strrep(bar_names, ['_' m{:}],'');
    end
    bar_names = unique(strrep(bar_names, '_', ' '));
    
    % remove sobt formulas from output
    bar_names(strcmp(bar_names,'sobt fv')) = [];
    bar_names(strcmp(bar_names,'sobt p')) = [];
    bar_names(strcmp(bar_names,'sobt pv')) = [];
    bar_names(strcmp(bar_names,'sobt so')) = [];
    bar_names(strcmp(bar_names,'sobt v')) = [];
    bar_names(strcmp(bar_names,'sobt vp')) = [];
    
    bar_cats = categorical(unique(bar_names));
    
    figure
    bar(bar_cats, cell2mat(ctime_mor))
    if strcmp(bench,'plate_48_rayleigh') || strcmp(bench,'plate_48_hysteretic')
        legend(["osimaginput", "osrealinput"]);
    else
        legend(proj_methods)
    end
    set(gca,'YScale','log')
    title(strrep(bench,'_',' '))
end

function this_names = typeset_names(this_names, proj_method)
        this_names = strrep(this_names, ['_' proj_method],'');
        this_names = strrep(this_names, '_', ' ');
        this_names = strrep(this_names, 'strint ', '');
        this_names = strrep(this_names, 'minrel', '\lmorminrel{}');
        this_names = strrep(this_names, 'avg', '\lmoravg{}');
        this_names = strrep(this_names, 'equi', '\lmorequi{}');
        this_names = strrep(this_names, 'linf', '\lmorlinf{}');
        this_names = strrep(this_names, 'strprs', '\lmorsp{}');
        this_names = strrep(this_names, 'aaaa', '\lmoraaa{}');
        this_names = strrep(this_names, 'std', '\lmorstd{}');
end

