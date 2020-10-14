function initialize_model(Info,MVars,StartDate,EndDate,SubDomain)
% Initialize required parameters for the output program
%
% Inputs: Info - Structure containing information about the model run
%         MVars - Structure containing the names of the model variables on
%           the output tape
%         StartDate - Starting Date for the simulation
%         EndDate - Ending Date for the simulation
%         SubDomain - Model Domain name
% No outputs, calls the routines to run the model
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated April 2017.

Info.StartDate = StartDate;
Info.EndDate = EndDate;
Info.NameIdentifier = SubDomain;
[modelpars] = model_pars;

% Right now, these options are turned off
if ~isfield(Info,'SaveForcingInterpFiles')
    Info.SaveForcingInterpFiles = 0;        % 0: Do not save daily forcing interpolation files
                                            % 1: Save daily forcing interpolation files (this takes
                                            % disk space but is much faster if doing repeated runs
                                            % for the same period)
end
if ~isfield(Info,'OnlyKeepLatestRestart')
    Info.OnlyKeepLatestRestart = 1;         % 0: Do not delete modeled restart files (keep all daily files)
                                         	% 1: Only keep the latest daily restart file
end
if ~isfield(Info,'ContinuationRun')
    Info.ContinuationRun = 0;               % 0: Normal model run
                                            % 1: Continuation run (If there is model output on a given
                                            % day, don't run for that day, otherwise attempt to load a
                                            % model state file and continue where the model left off
end
if ~isfield(Info,'InitializeFromRestartFile')
    Info.InitializeFromRestartFile = 1;    	% 0: Do not initialize a model run from an existing restart file
                                            % 1: Initialize a model run from an existing restart file (if it exists)
end

% Set the exponent for the IDW for the forcing correction interpolation (if used)
Info.forcingCorrExponent = 2;

% Constants
gen.MM2IN = 25.4;
gen.DAY = 86400;
gen.HOUR = 3600;
gen.KELVIN = 273.15;
gen.M2MM = 1000;
gen.TINY = 1E-6;
gen.LARGE = 1E6;
gen.TS = gen.HOUR;         % Model timestep (sec)
                                       
% Hard wired parameters (to make these user definable, move to
% model_pars.m) - Not intended for modification

% Constant energy terms
modelpars.fusheat       = 3.34e5;                  	% latent heat of fusion [J/kg]  
modelpars.subheat       = 2.85e6;                   % latent heat of sublimation [J/kg]
modelpars.evapheat      = 2.26e6;                   % latent heat of evaporation [J/kg]
modelpars.specheat_a    = 1008;                   	% specific heat of air [J/kg-K]
modelpars.specheat_w    = 4181;                     % specific heat of water [J/kg-K] 
modelpars.specheat_i    = 2050;                     % specific heat of ice [J/kg-K]
modelpars.specheat_s    = 1480;                     % specific heat of soil
modelpars.rhow          = 1000;                   	% density of water (kg/m3)
modelpars.rhoa          = 1.229;                  	% density of air   (kg/m3)
modelpars.rhoi          = 931;                      % density of ice   (kg/m3)
modelpars.rhos          = 1300;                     % density of soil   (kg/m3)
modelpars.g             = 9.81;                     % Gravity at earth's surface [m/s2]
modelpars.karman        = 0.41;                     % von karman constant
modelpars.emiss_snow    = 0.99;                     % Snow surface emmissivity
modelpars.boltz         = 5.6704e-8;              	% stefan-boltzmann constant
modelpars.kappa         = 0.514;                  	% soil thermal conductivity (Wm-1K-1)
modelpars.R_v           = 461.5;                    % water vapor gas constant [J/kg-K]
modelpars.emiss_tree    = 0.96;                     % emssivity of trees
modelpars.emiss_ground  = 0.94;                     % surface emssivity

% Vegetation Dynamics (for now, turn off), 
modelpars.maxLAI_c = modelpars.LAI_max;                         % Weekly Temperature below which ground vegetation is in its winter state [C]
modelpars.minLAI_c = modelpars.LAI_max;                         % Weekly Temperature above which ground vegetation is in its summer state [C]
modelpars.maxLAI_g = modelpars.LAI_max;                         % Minimum LAI (for dense ground vegetation) [-]
modelpars.minLAI_g = modelpars.LAI_max;                         % Maximum LAI (for dense ground vegetation) [-]
modelpars.LAI_c_tmax = 0;                       % Weekly Temperature below which canopy is in its winter state [C]
modelpars.LAI_c_tmin = 0;                       % Weekly Temperature above which canopy is in its summer state [C]
modelpars.LAI_g_tmax = 0;                       % Minimum LAI (for dense canopy) [-]
modelpars.LAI_g_tmin = 0;                       % Maximum LAI (for dense canopy) [-]

% Adjustments to forcing data (no adjustments for now)
modelpars.RainMult = modelpars.PrecipMult;   	% Rainfall multiplier [-] (make equal to snow multiplier)
modelpars.SnowMult = modelpars.PrecipMult;    	% Snowfall multiplier [-] - from Sagehn
modelpars.LRadMult = 1;                         % Longwave radiation multiplier [-]
modelpars.SRadMult = 1;                         % Shortwave radiation multiplier [-]

% Other multipliers (no adjustments for now)
modelpars.HMult = 1;
modelpars.EMult = 1;
modelpars.acsub_mult = 1;
modelpars.TempAdjustment = 0;

% Temperature and wind measurement level
modelpars.templevel = 2;                        % Height of temperature measurement [m]
modelpars.windlevel = 10;                       % Height of windspeed measurement [m]

% Snow Density and Snow Albedo
modelpars.d_par = modelpars.d_par/100;
modelpars.ripe_d_par = modelpars.ripe_d_par/100;
modelpars.maxage_d = (modelpars.maxdensity-modelpars.density_i)/modelpars.d_par;
modelpars.albedo_snow_reset = 1;                % Snowfall needed to reset the albedo
modelpars.ScatterEff = 0.00;                    % Re-irradiation of a snowpack due to scattering [-] - None for now

% Miscellaneous parameters (non snow)
modelpars.groundlayer_thickness = 0.1;          % Thickness of ground layer beneath snow
modelpars.initialSMC = 0.2;                     % Initial Soil Moisture
modelpars.albedo_ground = 0.2;                  % Canopy Albedo [-] - Shuttleworth, 2012
modelpars.kappa_d_tree = 100;                   % Kappa/d for trees (thermal conductivity/distance)
modelpars.svroughness = 0.05;                   % Roughness of bare ground surface [m]
modelpars.rc = 105;                             % Canopy Resistance [-]
modelpars.vzeroplan = 0.05;                     % Zero plane displacement height of vegetation elements above surface [m]
modelpars.vroughness = 0.1;                     % Roughenss of vegetateion elements above surface [m]

disp(['Running Model for "' Info.NameIdentifier]);
[Info,Cutout] = get_bounds(Info);

if Info.RunParallel == 1 && numel(Cutout) > 1 && Info.MaxNumWorkers > 1          % Built in Parallelization
    try 
        parpool(min(Info.MaxNumWorkers,numel(Cutout))); 
    end
    parfor i=1:numel(Cutout)
        Info2 = Info;
        Info2.NameIdentifier_orig = Info.NameIdentifier;
        Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
        Info2.Cutout = Cutout(i);
        IdentifierList(i).Identifier = Info2.NameIdentifier;
        [SpatialData] = process_spatial(Info2,modelpars);
        if ~isempty(SpatialData)
            [SpatialData] = POIs(Info2,SpatialData);
            SnowPALM_run(gen,Info2,SpatialData,modelpars,MVars);  
        end
    end
else                                % Serial Mode
    for i=1:numel(Cutout)
        Info2 = Info;
        Info2.NameIdentifier_orig = Info.NameIdentifier;
        Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
        Info2.Cutout = Cutout(i);
        IdentifierList(i).Identifier = Info2.NameIdentifier;
        [SpatialData] = process_spatial(Info2,modelpars);
        if ~isempty(SpatialData)
            [SpatialData] = POIs(Info2,SpatialData);
            SnowPALM_run(gen,Info2,SpatialData,modelpars,MVars);  
        end
    end
end
