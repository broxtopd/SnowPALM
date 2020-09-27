function modelpars = model_pars

% Adjustments to forcing dataset if necessary
modelpars.PrecipMult = 0.9;                       % Snowfall multiplier [-] - from Sagehn
modelpars.tsnow = 2;                            % Temperature at which precipitation is considered snow [C]

% Canopy Parameters
modelpars.LAI_max = 4;                          % LAI of fully dense forest [m2 m-2]
modelpars.snow_unload_par = 0.083;              % Snow unloading coefficient [day-1]
modelpars.albedo_tree = 0.15;                   % Canopy Albedo [-] - Shuttleworth, 2012
modelpars.EdgenessImportance = 0.7;            	% Describes how canopy radiation is transferred from the canopy to the ground                                       
                                                % Modified by EFLength in the horizontal and EFLength_vert in the vertical
modelpars.EFLength = 5;                         % Horizontal area of influence of vegetation skin [m] - integer (0-10)
modelpars.EFLength_vert = 10;                  	% Vertical area of influence of vegetation skin [m]

% Snowpack energy balance
modelpars.sroughness = 0.0003;           	    % Snow surface roughness length [m]
modelpars.kappa_s_density_i = 0.2;              % Thermal conductivity of new snow [W m-1 K-1]
modelpars.kappa_s_max_density = 0.7;            % Thermal conductivity of dense snow [W m-1 K-1]
modelpars.liq_swe_cap_mult = 0.033;             % Liquid water capacity of snow [-]

% Albedo and snow density
modelpars.albedo_i = 0.7;                     	% Albedo of fresh snow [-] (nadir looking)
modelpars.minalbedo = 0.3;                   	% Minimum possible snow albedo [-] (nadir looking)
modelpars.maxage_a = 10;                        % Albedo/Age Decay Constant [day] 
modelpars.ripe_a_par = 0.05;                 	% Albedo Decay (1/day) at at total liquid swe capacity [1/day]
modelpars.groundveght = 0.2;                    % Ground vegetation height [m] - estimate
modelpars.density_i = 0.1;                      % Density of fresh snow [mm swe/mm uncompacted snow]
modelpars.maxdensity = 0.55;                  	% Maximum snow density [mm swe/mm uncompacted snow]
modelpars.maxage_d = 70;                    	% Age that snow becomes fully old [day] 
modelpars.d_par = 0.65;                         % Decay constant for cold snow [%/day]
modelpars.ripe_d_par = 5;                       % Decay constant for warm snow [%/day]

% Wind Factors for wind redistribution metric (Winstrel et al.)
% Determine value of Sx index
modelpars.PrevailingWindAngle = 210;            % Prevailing Wind Angle (degrees from east, measured counterclockwise)
modelpars.Sx_Exponent = 10;                     % Exponent applied to Sx Index
modelpars.Sx_Rescale = 0;                   	% Minimum value of Sx (this value becomes zero)
modelpars.Sb_Sepdist = 100;                     % Separation distance for Sb
modelpars.Sb_Multiplier = 2;                    % Multiplier applied to Sb index

% Soil Parameters
modelpars.kappa_soil = 0.5;                     % Thermal conductivity of soil [W m-1 K-1]
modelpars.dampd = 2;                            % Damping depth (for ground heat flux below soil layer) [m]
modelpars.tempdampd = 2;                      	% Temperature at damping depth [C]

