function initialize_output(Info,FVars,MVars,StartDate,EndDate,SubDomain)
% Initialize required parameters for the output program
%
% Inputs: Info - Structure containing information about the model run
%         FVars - Structure containing the names of forcing variables on
%           the output tape
%         MVars - Structure containing the names of the model variables on
%           the output tape
%         StartDate - Starting Date for the simulation
%         EndDate - Ending Date for the simulation
%         SubDomain - Model Domain name
% No outputs, calls the routines to create output spatial data
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated April 2017.

Info.StartDate = StartDate;
Info.EndDate = EndDate;
[modelpars] = model_pars;

% Constants
gen.MM2IN = 25.4;
gen.DAY = 86400;
gen.HOUR = 3600;
gen.KELVIN = 273.15;
gen.M2MM = 1000;
gen.TINY = 1E-6;
gen.LARGE = 1E6;
gen.TS = gen.HOUR;         % Model timestep (sec)

% Don't alter the forcing
modelpars.RainMult = 1;                         % Rainfall multiplier [-]
modelpars.TempAdjustment = 0;                   % Temperature adjustment (C)
modelpars.LRadMult = 1;                         % Longwave radiation multiplier [-]
modelpars.SRadMult = 1;                         % Shortwave radiation multiplier [-]

modelpars.SnowMult = modelpars.PrecipMult;
modelpars.RainMult = modelpars.PrecipMult;

% Vegetation Dynamics (for now, turn off), 
modelpars.maxLAI_c = modelpars.LAI_max;                         % Weekly Temperature below which ground vegetation is in its winter state [C]
modelpars.minLAI_c = modelpars.LAI_max;                         % Weekly Temperature above which ground vegetation is in its summer state [C]
modelpars.maxLAI_g = modelpars.LAI_max;                         % Minimum LAI (for dense ground vegetation) [-]
modelpars.minLAI_g = modelpars.LAI_max;                         % Maximum LAI (for dense ground vegetation) [-]
modelpars.LAI_c_tmax = 0;                       % Weekly Temperature below which canopy is in its winter state [C]
modelpars.LAI_c_tmin = 0;                       % Weekly Temperature above which canopy is in its summer state [C]
modelpars.LAI_g_tmax = 0;                       % Minimum LAI (for dense canopy) [-]
modelpars.LAI_g_tmin = 0;                       % Maximum LAI (for dense canopy) [-]

Info.NameIdentifier = SubDomain;

% Output spatial data from each tile individually

% Get the tile bounds
[Info,Cutout] = get_bounds(Info);
CreateBoundsShapefile(Info,Cutout)

% Run the model (If doing more than one tile, and if the paralell options
% are on, use the built in paralellization)
% Built in Paralellization
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
        display_output(gen,Info2,modelpars,FVars,MVars);   
    end
    pause(2)
% Serial Mode 
else  
    for i=1:numel(Cutout)
        Info2 = Info;
        Info2.NameIdentifier_orig = Info.NameIdentifier;
        Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
        Info2.Cutout = Cutout(i);
        IdentifierList(i).Identifier = Info2.NameIdentifier;
        display_output(gen,Info2,modelpars,FVars,MVars); 
    end
    pause(2)
end

% Put the maps together
Info.Cutout = Cutout;
combine_maps(gen,Info,IdentifierList,FVars,MVars);
