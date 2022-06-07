function model = load_model(datfile, load_symfuns)
%LOAD_MODEL Loads the benchmark example DATFILE + stijns example

%
% This file is part of the Code, Data and Results for Numerical Experiments
% in "Structured model order reduction for vibro-acoustic problems using
% interpolation and balancing methods"
% Copyright (C) 2022 Quirin Aumann and Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin == 1
    load_symfuns = false;
end

% Location of model data
folder = 'models';

% Determine model type
if startsWith(datfile,'poroacoustic')
    type = 'poro_acoustic';
elseif startsWith(datfile,'radiation')
    type = 'acoustic';
elseif startsWith(datfile,'plate_48_hysteretic')
    type = 'hysteretic';
elseif startsWith(datfile,'plate_48_rayleigh') ...
        || startsWith(datfile,'transmission')
    type = 'proportional';
elseif startsWith(datfile,'test')
    type = 'proportional';
end

% Load data.
model = load([folder filesep() datfile '.mat']);
model.name = datfile;

% Align data
if isfield(model, 'B')
    model.b = model.B;
    model = rmfield(model,'B');
end

if isfield(model, 'C')
    model.c = model.C;
    model = rmfield(model,'C');
end

if isfield(model, 'c')
    if size(model.c, 1) > size(model.c, 2)
        model.c = transpose(model.c);
    end
end

% Arrange data
switch type
    case 'poro_acoustic'
        % Material parameters (from Rumpler 2014).
        density_fluid           = 1.21;
        viscosity_fluid         = 1.84e-5;
        standard_pressure_fluid = 101e3;
        heat_capacity_fluid     = 1.4;
        prandtl_number_fluid    = 0.71;
        porosity                = 0.96;
        tortuosity              = 1.7;
        flow_resistivity        = 32e3;
        viscous_length          = 90e-6;
        thermal_length          = 165e-6;

        % Frequency-dependent functions.
        viscous_drag = @(s) flow_resistivity * porosity^2 ...
            * sqrt(1 + 4 * s * tortuosity^2 * viscosity_fluid * density_fluid / ...
            (flow_resistivity^2 * viscous_length^2 * porosity^2));

        apparent_mass_density = porosity * density_fluid * (tortuosity - 1);

        alpha = @(s) 1 + 8 * viscosity_fluid / (s * prandtl_number_fluid ...
            * thermal_length^2 * density_fluid) ...
            * sqrt(1 + s * prandtl_number_fluid * thermal_length^2 ...
            * density_fluid / (16 * viscosity_fluid));

        fun{1} = @(s) 1;
        fun{2} = @(s) -porosity * (-apparent_mass_density - viscous_drag(s) / s) ...
            / (porosity * density_fluid + apparent_mass_density ...
            + viscous_drag(s) / s);
        fun{3} = @(s) porosity^2 / (porosity * density_fluid ...
            + apparent_mass_density + viscous_drag(s) / s);
        fun{4} = @(s) s^2;
        fun{5} = @(s) s^2 * (viscous_drag(s) / s - ...
            (-apparent_mass_density - viscous_drag(s) / s)^2 / ...
            (porosity * density_fluid + apparent_mass_density ...
            + viscous_drag(s) / s));
        fun{6} = @(s) s^2 * (-(heat_capacity_fluid-1) / (alpha(s) ...
            * heat_capacity_fluid * standard_pressure_fluid));
        fun{7} = @(s) s^2 * fun{2}(s);

        funB = @(s) -s^2;
        
        model.fA = fun;
        model.fb = funB;
        
        % load symbolic functions if required
        if load_symfuns
            
            syms vds(s) as(s)
            vds(s) = flow_resistivity * porosity^2 * sqrt(1 + 4*s*tortuosity^2*viscosity_fluid*density_fluid / ...
                (flow_resistivity^2*viscous_length^2*porosity^2));
            as(s) = 1+8*viscosity_fluid / (s*prandtl_number_fluid * ...
                thermal_length^2 * density_fluid) * sqrt(1+s*prandtl_number_fluid * ...
                thermal_length^2 * density_fluid / (16*viscosity_fluid));

            sfun{1}(s) = s/s;
            sfun{2}(s) = -porosity * (-apparent_mass_density - vds(s)/s) / ...
                (porosity*density_fluid + apparent_mass_density + vds(s)/s);
            sfun{3}(s) = porosity^2 / (porosity*density_fluid + apparent_mass_density + ...
                vds(s)/s);
            sfun{4}(s) = s^2;
            sfun{5}(s) = s^2 * (vds(s)/s - (-apparent_mass_density - vds(s)/s)^2 / ...
                (porosity*density_fluid + apparent_mass_density + vds(s)/s));
            sfun{6}(s) = s^2 * (-(heat_capacity_fluid-1) / (as(s) * heat_capacity_fluid*standard_pressure_fluid));
            sfun{7}(s) = s^2 * sfun{2}(s);

            sfunb(s) = -s^2;

        end
        
        if load_symfuns
            order = 6;
            % compute analytic derivatives
            model.symfuns.dsfuns = cell(order+1,length(sfun));
            for jj=1:length(sfun)
                model.symfuns.dsfuns{1,jj} = sfun{jj};
                for ii=1:order
                    model.symfuns.dsfuns{ii+1,jj}(s) = diff(sfun{jj}(s),s,ii);
                end
            end

            model.symfuns.dsfunb = cell(order+1,1);
            model.symfuns.dsfunb{1}(s) = sfunb(s);
            for ii=1:order
                model.symfuns.dsfunb{ii+1}(s) = diff(sfunb(s),s,ii) / factorial(ii);
            end
        end
        
    case 'acoustic'
        % we choose output 5 in this example
        n = 5;
        model.c = model.c(n,:);
        
        fun{1} = @(s) 1;
        fun{2} = @(s) s^2;
        funB = @(s) -s;
        
        model.fA = fun;
        model.fb = funB;
        
        % load symbolic functions if required
        if load_symfuns
            syms s
            
            model.symfuns.dsfuns = cell(3,2);
            model.symfuns.dsfuns{1,1}(s) = s/s;
            model.symfuns.dsfuns{2,1}(s) = 0*s;
            model.symfuns.dsfuns{3,1}(s) = 0*s;
            model.symfuns.dsfuns{1,2}(s) = s^2;
            model.symfuns.dsfuns{2,2}(s) = 2*s;
            model.symfuns.dsfuns{3,2}(s) = 1*(s/s);
            
            model.symfuns.dsfunb = cell(3,1);
            model.symfuns.dsfunb{1}(s) = -s;
            model.symfuns.dsfunb{2}(s) = -s/s;
            model.symfuns.dsfunb{3}(s) = 0*s;
        end
        
    case 'hysteretic'
        fun{1} = @(s) 1;
        fun{2} = @(s) s^2;
        
        model.fA = fun;
        model.fb = @(s) 1;
        
        % load symbolic functions if required
        if load_symfuns
            syms s
            
            model.symfuns.dsfuns = cell(3,2);
            model.symfuns.dsfuns{1,1}(s) = s/s;
            model.symfuns.dsfuns{2,1}(s) = 0*s;
            model.symfuns.dsfuns{3,1}(s) = 0*s;
            model.symfuns.dsfuns{1,2}(s) = s^2;
            model.symfuns.dsfuns{2,2}(s) = 2*s;
            model.symfuns.dsfuns{3,2}(s) = 1*(s/s);
        end
        
    case 'proportional'
        fun{1} = @(s) 1;
        fun{2} = @(s) s;
        fun{3} = @(s) s^2;
        
        model.fA = fun;
        model.fb = @(s) 1;
        
        % load symbolic functions if required
        if load_symfuns
            syms s
            
            model.symfuns.dsfuns = cell(3,3);
            model.symfuns.dsfuns{1,1}(s) = s/s;
            model.symfuns.dsfuns{2,1}(s) = 0*s;
            model.symfuns.dsfuns{3,1}(s) = 0*s;
            model.symfuns.dsfuns{1,2}(s) = s;
            model.symfuns.dsfuns{2,2}(s) = s/s;
            model.symfuns.dsfuns{3,2}(s) = 0*s;
            model.symfuns.dsfuns{1,3}(s) = s^2;
            model.symfuns.dsfuns{2,3}(s) = 2*s;
            model.symfuns.dsfuns{3,3}(s) = 1*(s/s);
            
        end
                
    otherwise
        error('unknown benchmark problem')
end

end
