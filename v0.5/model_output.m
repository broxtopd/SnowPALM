function model_output(StartDate,EndDate,SubDomain)
% Displays user requested maps or timeseries plots of model output or
% spatial data

%% Begin User Input
Info.RunParallel = 1;               % 0: No Paralellization (run each chunk sequentially), 1: Use Matlab's built in Paralellzation
Info.MaxNumWorkers = 16;            % Maximum number of workers to work on modelling simultaneously

% Basic Terrain Informaion to output
Info.DEM         = 1;           % Digital Elavation Model    
Info.Slope       = 1;           % Slope Map
Info.Northness   = 1;           % Northness Map   
Info.VegHT       = 1;           % Vegetation Height Map
Info.Cover       = 1;           % Overhead Canopy Cover Map
Info.SkyView     = 1;           % Sky View Factor Map

% Composite Model Indicies to output
Info.SFI_BareEarth   = 1;    	% SFI Bare Earth
Info.SFI_Veg         = 1;    	% SFI (with vegetation)
Info.Sx         	 = 1;       % Sx_Index (-)
Info.Sb              = 1;       % Sb_Index (-)

Info.outputType = 'day';            % Output on daily or hourly timestep (hour, day)
Info.outputFType = 'tif';           % Output file type (tif, mat) 
Info.CreateOverviewTiles = 1;       % Create pyramids for GIS output

FVars(1).Name = 'Rain';             FVars(1).map = 0;        % Rainfall (mm)   (compute from forcing data)
FVars(2).Name = 'Snow';             FVars(2).map = 0;        % Snowfall (mm)   (compute from forcing data)
FVars(3).Name = 'AirT';             FVars(3).map = 0;        % Air Temperature (K)   (compute from forcing data)
FVars(4).Name = 'Pres';             FVars(4).map = 0;        % Pressure (Pa)   (compute from forcing data)
FVars(5).Name = 'RH';               FVars(5).map = 0;        % Relative Humidity (%)   (compute from forcing data)
FVars(6).Name = 'TDew';             FVars(6).map = 0;        % Dew Point Temperature (C)   (compute from forcing data)
FVars(7).Name = 'Vapp';             FVars(7).map = 0;     	 % Vapor Pressure (Pa)   (compute from forcing data)
FVars(8).Name = 'Wind';             FVars(8).map = 0;        % Wind Speed (m/s)   (compute from forcing data)
FVars(9).Name = 'WindAngle';        FVars(9).map = 0;        % Wind Angle (degrees from east, measured counterclockwise)   (compute from forcing data)
FVars(10).Name = 'DSWRF';           FVars(10).map = 0;       % Downward Shortwave Radiation (W/m2)   (compute from forcing data)
FVars(11).Name = 'DLWRF';           FVars(11).map = 0;       % Downward Longwave Radiation (W/m2)   (compute from forcing data)

MVars(1).Name = 'Rain';                    MVars(1).map = 0;          % Rainfall (mm/hr) 
MVars(2).Name = 'Snow';                    MVars(2).map = 0;          % Snowfall (mm/hr)
MVars(3).Name = 'AirT';                    MVars(3).map = 0;          % Air Temperature (C)
MVars(4).Name = 'Pres';                    MVars(4).map = 0;          % Pressure (Pa)
MVars(5).Name = 'RH';                      MVars(5).map = 0;          % Relative Humidity (%)
MVars(6).Name = 'TDew';                    MVars(6).map = 0;          % Dew Point Temperature (C)
MVars(7).Name = 'Vapp';                    MVars(7).map = 0;          % Vapor Pressure (Pa)
MVars(8).Name = 'Wind';                    MVars(8).map = 0;          % Wind Speed (m/s)
MVars(9).Name = 'WindAngle';               MVars(9).map = 0;          % Wind Angle (degrees from east, measured counterclockwise)
MVars(10).Name = 'DSWRF';                  MVars(10).map = 0;         % Downward Shortwave Radiation (W/m2)
MVars(11).Name = 'DLWRF';                  MVars(11).map = 0;         % Downward Longwave Radiation (W/m2)
MVars(12).Name = 'Adj_Snowfall';           MVars(12).map = 0;         % Snowfall (After accounting for undercatch and wind effects)       (mm/hr)    
MVars(13).Name = 'TFall_Snow';             MVars(13).map = 0;         % Throughfall (snow - mm/hr)      
MVars(14).Name = 'TFall_Rain';             MVars(14).map = 0;         % Throughfall (rain - mm/hr)
MVars(15).Name = 'Snow_Unload';            MVars(15).map = 0;         % Canopy Snow unloading (mm/hr)      
MVars(16).Name = 'Melt_Drip';              MVars(16).map = 0;         % Canpy melt drip (mm/hr)  
MVars(17).Name = 'Canopy_Sublimation';     MVars(17).map = 0;         % Canopy sublimation (mm/hr)   
MVars(18).Name = 'Canopy_Snow_Storage';    MVars(18).map = 0;         % Canopy snow storage (mm)       
MVars(19).Name = 'SWE';                    MVars(19).map = 0;         % SWE (mm)       
MVars(20).Name = 'Snowpack_Liquid';        MVars(20).map = 0;         % Liquid in the snowpack (mm)     
MVars(21).Name = 'Snowpack_Sublimation';   MVars(21).map = 0;         % Snowpack Sublimation (mm/hr)     
MVars(22).Name = 'Snow_Melt';              MVars(22).map = 0;         % Snow Melt (mm/hr)     
MVars(23).Name = 'Snow_Density';           MVars(23).map = 0;         % Snow Density (m3/m3)
MVars(24).Name = 'Snow_Depth';             MVars(24).map = 0;         % Snow Depth (mm)
MVars(25).Name = 'Snow_Temp_s';            MVars(25).map = 0;         % Snow surface layer temperature (C)      
MVars(26).Name = 'Snow_Temp_m';            MVars(26).map = 0;         % Temperature of bulk snowpack (C) 
MVars(27).Name = 'Albedo_Snow';            MVars(27).map = 0;         % Snow Albedo (-)
MVars(28).Name = 'Net_Radiation_Snow';     MVars(28).map = 0;         % Net Radiation at snow surface (W/m2)   
MVars(29).Name = 'Sensible_Heat_Snow';     MVars(29).map = 0;         % Sensible heat flux from the snow surface (W/m2)    
MVars(30).Name = 'Latent_Heat_Snow';       MVars(30).map = 0;         % Latent heat flux from the snow surface (W/m2)    
MVars(31).Name = 'Heat_From_Precip_Snow';  MVars(31).map = 0;       	% Heat flux from falling precip (W/m2)   
MVars(32).Name = 'Ground_Heat';            MVars(32).map = 0;       	% Ground heat flux (W/m2)   
MVars(33).Name = 'Cold_Content';           MVars(33).map = 0;         % Snowpack cold content (J/m2)     
MVars(34).Name = 'Residual_To_Soil';       MVars(34).map = 0;         % Energy imbalence term (goes to soil layer) (W/m2)
MVars(35).Name = 'Heat_From_Phase_Change'; MVars(35).map = 0;         % Heat associated with melting and refreezing of snowpack (W/m2)        
MVars(36).Name = 'Net_Radiation_Ground';   MVars(36).map = 0;         % Net Radiadiation at ground surface (W/m2)
MVars(37).Name = 'SRad_in_Ground';         MVars(37).map = 0;         % Incoming longwave radiation to ground surface (W/m2)
MVars(38).Name = 'LRad_in_Ground';         MVars(38).map = 0;         % Incoming shortwave radiation to ground suraface (W/m2) (W/m2)
MVars(39).Name = 'SRad_Scatter_to_Ground'; MVars(39).map = 0;         % Shortwave scattering to ground surface (W/m2)         
MVars(40).Name = 'Rc_add_g';               MVars(40).map = 0;         % Energy imbalance added to ground surface (W/m2)  
MVars(41).Name = 'Ground_Albedo';          MVars(41).map = 0;         % Ground Albedo        
MVars(42).Name = 'Ground_Emissivity';      MVars(42).map = 0;         % Ground Emmissivity
MVars(43).Name = 'Temp_Ground';            MVars(43).map = 0;         % Ground surface temperature (C)                 
MVars(44).Name = 'LAI';                    MVars(44).map = 0;         % Leaf area index (m3/m2)
MVars(45).Name = 'Q_soil';                 MVars(45).map = 0;         % Energy content of soil layer (J/m2)       
MVars(46).Name = 'ice_soil';               MVars(46).map = 0;         % Ice content of soil layer (mm)      
MVars(47).Name = 'T_soil';                 MVars(47).map = 0;         % Temperature of soil layer (C) 
MVars(48).Name = 'EFlux_soil';             MVars(48).map = 0;         % Soil energy flux (W/m2)

%% End User Input

[Info] = program_pars(SubDomain,Info);   % Get spatial parameters
PWD = pwd;
cd(Info.ProgramFilesDir)
set_paths(Info.ProgramFilesDir,'add')
cd(PWD)

% Only compute the 3-D indices if the wind or sky view maps are to be made
if ~strcmpi(Info.SxSource,'None')
    if (Info.Sx || Info.Sb) && ~strcmpi(Info.SxSource,'None')
        Info.ComputeSxIndex = 1;
    else
        Info.ComputeSxIndex = 0;
    end
    if (Info.Sb) && ~strcmpi(Info.SxSource,'None')
        Info.ComputeSbIndex = 1;
    else
        Info.ComputeSbIndex = 0;
    end
else
    Info.ComputeSxIndex = 0;
    Info.ComputeSbIndex = 0;
    if (Info.Sx || Info.Sb)
        disp('Warning: Wind indices won''t be computed because the SxSource variable is set to ''None''');
    end
end
if Info.SkyView || Info.Sx || Info.Sb || Info.SFI_BareEarth || Info.SFI_Veg
    Info.ComputeAngles = 1;
else
    Info.ComputeAngles = 0;
end

preprocess(Info,SubDomain)
initialize_output(Info,FVars,MVars,StartDate,EndDate,SubDomain)