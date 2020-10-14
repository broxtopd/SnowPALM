function run_model(StartDate,EndDate,SubDomain)
% Run the Snow Model
% Created by Patrick Broxton
% Updated 11/10/2012

% run_model('10/1/2009','11/1/2009','MConTower')

%% Begin User Input
      
Info.UseWindModel = 1;          % 1: Simulate Wind Redistribution, 0: Do not simulate wind redistribution 

Info.RunParallel = 1;           % 0: No Paralellization, 1: Built in Paralellzation
Info.MaxNumWorkers = 16;        % Maximum Number of paralell workers

Info.MapType = 'day';           % Flag to tell whether to output daily maps of maps every timestep (ts, day)
Info.SaveDailyRestartFiles = 0;
Info.SaveForcingInterpFiles = 0;    % Saves forcing interpolation files (much faster 

MVars(1).Name = 'Rain';                    MVars(1).modelvar = 'rain';                  MVars(1).map = 1;          % Rainfall (mm/hr) 
MVars(2).Name = 'Snow';                    MVars(2).modelvar = 'snow_0';             	MVars(2).map = 1;          % Snowfall (mm/hr)
MVars(3).Name = 'AirT';                    MVars(3).modelvar = 'airt';                  MVars(3).map = 1;          % Air Temperature (C)
MVars(4).Name = 'Pres';                    MVars(4).modelvar = 'pres';                  MVars(4).map = 1;          % Pressure (Pa)
MVars(5).Name = 'RH';                      MVars(5).modelvar = 'rh';                    MVars(5).map = 1;          % Relative Humidity (%)
MVars(6).Name = 'TDew';                    MVars(6).modelvar = 'tdew';                  MVars(6).map = 1;          % Dew Point Temperature (C)
MVars(7).Name = 'Vapp';                    MVars(7).modelvar = 'vapp';                  MVars(7).map = 1;          % Vapor Pressure (Pa)
MVars(8).Name = 'Wind';                    MVars(8).modelvar = 'wind';                  MVars(8).map = 1;          % Wind Speed (m/s)
MVars(9).Name = 'WindAngle';               MVars(9).modelvar = 'windangle';             MVars(9).map = 1;          % Wind Angle (degrees from east, measured counterclockwise)
MVars(10).Name = 'DSWRF';                  MVars(10).modelvar = 'dswrf';                MVars(10).map = 1;         % Downward Shortwave Radiation (W/m2)
MVars(11).Name = 'DLWRF';                  MVars(11).modelvar = 'dlwrf';                MVars(11).map = 1;         % Downward Longwave Radiation (W/m2)
MVars(12).Name = 'Adj_Snowfall';           MVars(12).modelvar = 'snow';                 MVars(12).map = 1;         % Snowfall (After accounting for undercatch and wind effects)       (mm/hr)
MVars(13).Name = 'TFall_Snow';             MVars(13).modelvar = 'tsfall';               MVars(13).map = 1;         % Throughfall (snow - mm/hr) 
MVars(14).Name = 'TFall_Rain';             MVars(14).modelvar = 'tfall';                MVars(14).map = 1;         % Throughfall (rain - mm/hr)
MVars(15).Name = 'Snow_Unload';            MVars(15).modelvar = 'snow_unload';          MVars(15).map = 1;         % Canopy Snow unloading (mm/hr) 
MVars(16).Name = 'Melt_Drip';              MVars(16).modelvar = 'melt_drip';            MVars(16).map = 1;         % Canpy melt drip (mm/hr)  
MVars(17).Name = 'Canopy_Sublimation';     MVars(17).modelvar = 'acsub';                MVars(17).map = 1;         % Canopy sublimation (mm/hr)  
MVars(18).Name = 'Canopy_Snow_Storage';    MVars(18).modelvar = 'state.cansnowstor';    MVars(18).map = 1;         % Canopy snow storage (mm)     
MVars(19).Name = 'SWE';                    MVars(19).modelvar = 'state.swe';            MVars(19).map = 1;         % SWE (mm)        
MVars(20).Name = 'Snowpack_Liquid';        MVars(20).modelvar = 'state.liq_swe';        MVars(20).map = 1;         % Liquid in the snowpack (mm)      
MVars(21).Name = 'Snowpack_Sublimation';   MVars(21).modelvar = 'sublimation';          MVars(21).map = 1;         % Snowpack Sublimation (mm/hr)   
MVars(22).Name = 'Snow_Melt';              MVars(22).modelvar = 'melt';                 MVars(22).map = 1;         % Snow Melt (mm/hr) 
MVars(23).Name = 'Snow_Density';           MVars(23).modelvar = 'snowdensity';          MVars(23).map = 1;         % Snow Density (m3/m3) 
MVars(24).Name = 'Snow_Depth';             MVars(24).modelvar = 'state.sdepth';         MVars(24).map = 1;         % Snow Depth (mm) 
MVars(25).Name = 'Snow_Temp_s';            MVars(25).modelvar = 'state.T_snows';        MVars(25).map = 1;         % Snow surface layer temperature (C)       
MVars(26).Name = 'Snow_Temp_m';            MVars(26).modelvar = 'state.T_snowm';        MVars(26).map = 1;         % Temperature of bulk snowpack (C)  
MVars(27).Name = 'Albedo_Snow';            MVars(27).modelvar = 'albedo';               MVars(27).map = 1;         % Snow Albedo (-) 
MVars(28).Name = 'Net_Radiation_Snow';     MVars(28).modelvar = 'Rn_snow';              MVars(28).map = 1;         % Net Radiation at snow surface (W/m2)    
MVars(29).Name = 'Sensible_Heat_Snow';     MVars(29).modelvar = 'H';                    MVars(29).map = 1;         % Sensible heat flux from the snow surface (W/m2)     
MVars(30).Name = 'Latent_Heat_Snow';       MVars(30).modelvar = 'lambdaE';              MVars(30).map = 1;         % Latent heat flux from the snow surface (W/m2)    
MVars(31).Name = 'Heat_From_Precip_Snow';  MVars(31).modelvar = 'A_p';                  MVars(31).map = 1;         % Heat flux from falling precip (W/m2)       
MVars(32).Name = 'Ground_Heat';            MVars(32).modelvar = 'g';                    MVars(32).map = 1;         % Ground heat flux (W/m2)    
MVars(33).Name = 'Cold_Content';           MVars(33).modelvar = 'state.Q_snow';         MVars(33).map = 1;         % Snowpack cold content (J/m2)    
MVars(34).Name = 'Residual_To_Soil';       MVars(34).modelvar = 'Q_soil_add';           MVars(34).map = 1;         % Energy imbalence term (goes to soil layer) (W/m2) 
MVars(35).Name = 'Heat_From_Phase_Change'; MVars(35).modelvar = 'Q_conv';               MVars(35).map = 1;         % Heat associated with melting and refreezing of snowpack (W/m2)            
MVars(36).Name = 'Net_Radiation_Ground';   MVars(36).modelvar = 'Rn_g';                 MVars(36).map = 1;         % Net Radiadiation at ground surface (W/m2) 
MVars(37).Name = 'SRad_in_Ground';         MVars(37).modelvar = 'Rs_down_g';            MVars(37).map = 1;         % Incoming longwave radiation to ground surface (W/m2) 
MVars(38).Name = 'LRad_in_Ground';         MVars(38).modelvar = 'Rl_g';                 MVars(38).map = 1;         % Incoming shortwave radiation to ground suraface (W/m2) (W/m2) 
MVars(39).Name = 'SRad_Scatter_to_Ground'; MVars(39).modelvar = 'Rs_scatter';           MVars(39).map = 1;         % Shortwave scattering to ground surface (W/m2)               
MVars(40).Name = 'Rc_add_g';               MVars(40).modelvar = 'Rc_add_g';             MVars(40).map = 1;         % Energy imbalance added to ground surface (W/m2)   
MVars(41).Name = 'Ground_Albedo';          MVars(41).modelvar = 'albedo_g';             MVars(41).map = 1;         % Ground Albedo    
MVars(42).Name = 'Ground_Emissivity';      MVars(42).modelvar = 'emiss_g';              MVars(42).map = 1;         % Ground Emmissivity 
MVars(43).Name = 'Temp_Ground';            MVars(43).modelvar = 'Temp_g';               MVars(43).map = 1;         % Ground surface temperature (C)
MVars(44).Name = 'LAI';                    MVars(44).modelvar = 'LAI_c';                MVars(44).map = 1;         % Leaf area index (m3/m2) 
MVars(45).Name = 'Q_soil';                 MVars(45).modelvar = 'Q_soil';               MVars(45).map = 1;         % Energy content of soil layer (J/m2)     
MVars(46).Name = 'ice_soil';               MVars(46).modelvar = 'ice_soil';             MVars(46).map = 1;         % Ice content of soil layer (mm)  
MVars(47).Name = 'T_soil';                 MVars(47).modelvar = 'T_soil';               MVars(47).map = 1;         % Temperature of soil layer (C)  
MVars(48).Name = 'EFlux_soil';             MVars(48).modelvar = 'p_dQ';                 MVars(48).map = 1;         % Soil energy flux (W/m2) 

%% End User Input

[Info] = program_pars(SubDomain,Info);   % Get spatial parameters
PWD = pwd;
cd(Info.ProgramFilesDir)
set_paths(Info.ProgramFilesDir,'add')
cd(PWD)

if Info.UseWindModel
    if ~strcmpi(Info.SxSource,'None')
        Info.ComputeSxIndex = 1;
        Info.ComputeSbIndex = 1;
    else
        disp('Warning: Wind model won''t be run because the SxSource variable is set to ''None''');
        Info.ComputeSxIndex = 0;
        Info.ComputeSbIndex = 0;
    end
else
    Info.ComputeSxIndex = 0;
    Info.ComputeSbIndex = 0;
end
Info.ComputeAngles = 1;         % For now, always use the horizon angle stuff

preprocess(Info,SubDomain)
initialize_model(Info,MVars,StartDate,EndDate,SubDomain)
output_csv(Info,StartDate,EndDate,SubDomain,MVars)