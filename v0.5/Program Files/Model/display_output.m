function display_output(gen,Info,modelpars,FVars,MVars)
% Print model outputs to GIS Maps
%
% Inputs: gen - Structure containing constants (e.g. seconds in model timestep)
%         Info - Structure containing information about the model run
%         modelpars - Structure containing model parameters
%         StartDate - Starting Date for the simulation
%         FVars - Structure containing the names of forcing variables on
%           the output tape
%         MVars - Structure containing the names of the model variables on
%           the output tape
% No outputs, calls the routines to create output spatial data
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)

%% Set Up Options
STStamp = datenum(Info.StartDate);
ETStamp = datenum(Info.EndDate);

% Get Spatial Data
[SpatialData] = process_spatial(Info,modelpars);

% Only operate on tiles that are defined, otherwise do nothing
if ~isempty(SpatialData)
    Z = SpatialData.z;
    S = SpatialData.S;
    R = SpatialData.R;
    veght = SpatialData.canopy_height;
    cover = SpatialData.canopy_coverage;

    % Only load forcing data if user requests a forcing map
    LoadForcingData = 0;
    for v = 1:numel(FVars)
        if FVars(v).map == 1
            LoadForcingData = 1;
        end
    end

    % If so, then compute variables needed for focing extraction program
    if LoadForcingData
        % Allow multiple forcing data sources
        for i = 1:numel(Info.ForcingDirs)
            % DEM, latitude, longitude of entire forcing grid (e.g. entire US)
            lon = ncread([Info.ForcingDirs{i} filesep 'dem.nc'],'lon');
            lat = ncread([Info.ForcingDirs{i} filesep  'dem.nc'],'lat');
            z = ncread([Info.ForcingDirs{i} filesep 'dem.nc'],'Elev');
            % Create 2d maps from vectors
            lon2d = repmat(lon,[1 numel(lat)]);
            lat2d = repmat(lat',[numel(lon) 1]);
            % Find which grids to extract
            xlocs = (lon >= Info.Forcing_ULLR(1)) & (lon <= Info.Forcing_ULLR(3));
            ylocs = (lat >= Info.Forcing_ULLR(4)) & (lat <= Info.Forcing_ULLR(2));
            % Latitude, Longitude, elevation grid
            Info.Forcing.lon2d{i} = lon2d(xlocs,ylocs);
            Info.Forcing.lat2d{i} = lat2d(xlocs,ylocs);
            Info.Forcing.Z{i} = z(xlocs,ylocs);
            % Size of extracted area
            Info.Forcing.sz{i} = size(z(xlocs,ylocs));
            % Positions in the forcing file where to extract data from
            Info.Forcing.start{i} = [find(lon >= Info.Forcing_ULLR(1),1,'first') find(lat >= Info.Forcing_ULLR(4),1,'first')];
            Info.Forcing.count{i} = [find(lon <= Info.Forcing_ULLR(3),1,'last')-Info.Forcing.start{i}(1)+1 find(lat <= Info.Forcing_ULLR(2),1,'last')-Info.Forcing.start{i}(2)+1];
            Info.Forcing.stride{i} = [1 1];
        end
        % Used when constructing file name to read from
        for j = 1:numel(Info.Forcing.var)
            ForcingDir = Info.ForcingDirs{Info.Forcing.source{j}};
            fprefix = strsplit(ForcingDir,filesep);
            fprefix = strsplit(char(fprefix(end)),'/');
            Info.Forcing.fprefix{j} = char(fprefix(end));
        end
    end

    %% General Terrain Maps
    fpath = [Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier];

    % These maps are based on data already in the SpatialData structure
    % (nothing to compute here)
    
    % Digital Elevation Model
    if Info.DEM   
        disp(['Getting DEM data for ' Info.NameIdentifier]);
        val = Z; 
        % Name of output map
        fname = 'elev';
        print_maps(Info,val,Z+veght,fpath,fname);
    end

    % Slope Map
    if Info.Slope
        disp(['Getting Slope data for ' Info.NameIdentifier]);
        val = S * 100; val(isnan(S)) = NaN;	
        fname = 'slope';
        print_maps(Info,val,Z+veght,fpath,fname);
    end

    % Northness Map
    if Info.Northness   
        disp(['Getting Northness data for ' Info.NameIdentifier]);
        val = sin(S) .* cos(R); 
        fname = 'northness';
        print_maps(Info,val,Z+veght,fpath,fname);
    end

    % Vegetation Height Map
    if Info.VegHT
        disp(['Getting Vegetation Height data for ' Info.NameIdentifier]);
        val = veght;
        fname = 'veght';
        print_maps(Info,val,Z+veght,fpath,fname);
    end

    % Vegetation Cover Map
    if Info.Cover 
        disp(['Getting Cover data for ' Info.NameIdentifier]);
        val = cover; 
        fname = 'cdensity';
        print_maps(Info,val,Z+veght,fpath,fname);
    end

    % Sky View
    if Info.SkyView
        vegcover_canopy = SpatialData.canopy_coverage/100;
        SkyView_noveg = 0;
        SkyView_canopy = 0;
        % Compute skyview as average skyview from each increment in
        % Alphamat (horizon angle in a partucular direction)
        for i = 1:numel(SpatialData.Alphas)-1
            SkyView_noveg = SkyView_noveg + (pi/2 - (double(SpatialData.AlphaMat(:,:,i))/(255/(pi/2)))) / (pi/2) / (numel(SpatialData.Alphas)-1); 
            SkyView_canopy = SkyView_canopy + (pi/2 - (double(SpatialData.VegAlphaMat(:,:,i))/(255/(pi/2)))) / (pi/2) / (numel(SpatialData.Alphas)-1); 
        end
        % For now, under canopy, scale from no coverage at zero canopy
        % height to complete coverage at 10 meter canopy height)
        adj_cover = max(max(0,min(1,veght/15)),vegcover_canopy);
        SkyView = min(SkyView_canopy,SkyView_noveg) .* (1-adj_cover).^2;

        val = SkyView; 
        fname = 'skyview';
        print_maps(Info,val,Z+veght,fpath,fname);
    end
        

    % Get wind index values if specified
    if Info.Sx || Info.Sb
        disp(['Computing Sx and Sb wind indexes for ' Info.NameIdentifier]);
        if Info.Sx
            sx = SpatialData.Sx;
            val = sx; 	
            disp(['Getting Sx data for ' Info.NameIdentifier]);
            fname = 'Sx';
            print_maps(Info,val,Z+veght,fpath,fname);
        end

        if Info.Sb
            sb = SpatialData.Sb * modelpars.Sb_Multiplier + 1;
            val = sb; 	
            disp(['Getting Sb data for ' Info.NameIdentifier]);
            fname = 'Sb';
            print_maps(Info,val,Z+veght,fpath,fname);
        end
    end

    % Solar forcing index based on Bare Earth Only
    if Info.SFI_BareEarth
        disp(['Computing Bare Earth Solar Forcing for ' Info.NameIdentifier]);
        [sfi] = getsolarindexes(Info,SpatialData,modelpars,0);
        val = sfi; 	
        fname = [datestr(datenum(Info.StartDate),'yyyymmdd') '_' datestr(datenum(Info.EndDate),'yyyymmdd')];
        fpath = [Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep 'sfi'];
        print_maps(Info,val,Z+veght,fpath,fname);
    end
    
    % Solar forcing index based on Bare Earth + Veg
    if Info.SFI_Veg
        disp(['Computing Solar Forcing (with canopy) for ' Info.NameIdentifier]);
        [sfi] = getsolarindexes(Info,SpatialData,modelpars,1);
        val = sfi; 	
        fname = [datestr(datenum(Info.StartDate),'yyyymmdd') '_' datestr(datenum(Info.EndDate),'yyyymmdd')];
        fpath = [Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep 'sfi_veg'];
        print_maps(Info,val,Z+veght,fpath,fname);
    end

    %% Model Variables on Particular Dates

    if ~(strcmp(Info.outputType,'avg') || strcmp(Info.outputType,'max') || strcmp(Info.outputType,'sum'))
        d = 0;
        for TS = STStamp:ETStamp
            % Get Forcing Maps (Day and TS)
            d = d+1;
            if LoadForcingData
                % Load Daily Forcing Data (this will either return daily
                % maps or print hourly maps)
                [rain_t,snow_t,airt_t,pres_t,rh_t,tdew_t,vapp_t,ugrd_t,vgrd_t,dswrf_t,dlwrf_t] = get_daily_forcing_data(gen,Info,FVars,modelpars,SpatialData,Z,veght,TS);
                % Compute a wind angle and total wind from the u and v wind
                % components
                windangle = atan2(vgrd_t,ugrd_t);
                windangle = windangle + pi;
                windangle_p = windangle + 2*pi;
                windangle_p(windangle > 0) = 0;
                windangle(windangle < 0) = 0;
                windangle_t = windangle + windangle_p;
                wind_t = sqrt(ugrd_t.^2+vgrd_t.^2);
                
                % For daily maps: Loop through forcing vars, check if they 
                % are requested, and make a map of them if they are
                for v = 1:numel(FVars)
                    if FVars(v).map && strcmp(Info.outputType,'day')
                        eval(['val = ' lower(FVars(v).Name) '_t;']);  
                        disp(['Getting ' FVars(v).Name ' data for ' datestr(TS,'yyyy-mm-dd') ' - ' Info.NameIdentifier]);
                        eval(['fname = datestr(TS,''yyyymmdd'');']);
                        eval(['fpath = [Info.DisplayDir filesep ''GIS'' filesep Info.NameIdentifier filesep ''' FVars(v).Name '_day''];']);
                        print_maps(Info,val,Z+veght,fpath,fname);
                    end
                end
            end

            % Get Model Maps (Day and TS)
            for v = 1:numel(MVars)
                if MVars(v).map
                    % Get daily model outputs if requested
                    if strcmp(Info.outputType,'day')
                        dirname = [Info.ModelStatesDir filesep 'LSM States' filesep 'Daily' filesep Info.NameIdentifier];
                        eval(['load([dirname filesep datestr(TS,''yyyymmdd'')],''' MVars(v).Name '_output'');']);
                        eval(['val = ' MVars(v).Name '_output;']); 
                        disp(['Getting ' MVars(v).Name ' data for ' datestr(TS,'yyyy-mm-dd') ' - ' Info.NameIdentifier]);
                        eval(['fname = datestr(TS,''yyyymmdd'');']);
                        eval(['fpath = [Info.DisplayDir filesep ''GIS'' filesep Info.NameIdentifier filesep ''' MVars(v).Name '_day''];']);
                        print_maps(Info,val,Z+veght,fpath,fname);
                    % Get hourly model outputs if requested
                    elseif strcmp(Info.outputType,'ts')
                        dirname = [Info.ModelStatesDir filesep 'LSM States' filesep 'TS' filesep Info.NameIdentifier];
                        eval(['load([dirname filesep datestr(TS,''yyyymmdd'')],''' MVars(v).Name '_output'');']);
                        for j = 1:gen.DAY/gen.TS
                            TS_sub = TS+((j-1)/(gen.DAY/gen.TS));
                            disp(['Getting ' MVars(v).Name ' data for ' datestr(TS_sub,'yyyy-mm-dd HH:MM') ' - ' Info.NameIdentifier]);
                            eval(['val = ' MVars(v).Name '_output(:,:,j);']);
                            eval(['fname = datestr(TS_sub,''yyyymmddHHMM'');']);
                            eval(['fpath = [Info.DisplayDir filesep ''GIS'' filesep Info.NameIdentifier filesep ''' MVars(v).Name '_ts''];']);
                            print_maps(Info,val,Z+veght,fpath,fname);
                        end
                    end
                end
            end
        end
    else
        %% Get Forcing Vars (Cumulative) - Not currently implemented, but theoretically works
        if LoadForcingData
            for v = 1:numel(FVars)
                eval([FVars(v).Name '_map = 0;']);
            end
            d = 0;
            STStamp = datenum(Info.StartDate);
            ETStamp = datenum(Info.EndDate);
            for TS = STStamp:ETStamp
                d = d+1;
                [rain_t,snow_t,airt_t,pres_t,rh_t,tdew_t,vapp_t,ugrd_t,vgrd_t,dswrf_t,dlwrf_t] = get_daily_forcing_data(gen,Info,FVars,modelpars,SpatialData,Z,veght,TS);
                windangle = atan2(vgrd_t,ugrd_t);
                windangle = windangle + pi;
                windangle_p = windangle + 2*pi;
                windangle_p(windangle > 0) = 0;
                windangle(windangle < 0) = 0;
                windangle_t = windangle + windangle_p;
                wind_t = sqrt(ugrd_t.^2+vgrd_t.^2);
                for v = 1:numel(FVars)
                    if FVars(v).map
                        disp(['Getting ' FVars(v).Name ' data for ' datestr(TS,'yyyy-mm-dd') ' - ' Info.NameIdentifier]);
                        eval([lower(FVars(v).Name) '_t(isnan(' lower(FVars(v).Name) '_t)) = 0;']);
                        if strcmp(Info.outputType,'avg') || strcmp(Info.outputType,'sum')
                            eval([FVars(v).Name '_map = ' FVars(v).Name '_map + ' lower(FVars(v).Name) '_t;']);
                        elseif strcmp(Info.outputType,'max')
                            eval([FVars(v).Name '_map = max(' FVars(v).Name '_map, ' lower(FVars(v).Name) '_t);']);
                        end
                    end
                end
            end

            for v = 1:numel(FVars)
                if FVars(v).map
                    if strcmp(Info.outputType,'avg')
                        eval([FVars(v).Name '_map = ' FVars(v).Name '_map / d;']);
                    end
                    eval(['val = ' FVars(v).Name '_map;']); 
                    fname = [Info.outputType '_' datestr(STStamp,'yyyymmdd') '_' datestr(ETStamp,'yyyymmdd')];
                    eval(['fpath = [Info.DisplayDir filesep ''GIS'' filesep Info.NameIdentifier filesep ''' FVars(v).Name '_composite''];']);
                    print_maps(Info,val,Z+veght,fpath,fname);
                end
            end
        end

        % Get Model Vars (Cumulative)
        for v = 1:numel(MVars)
            eval([MVars(v).Name '_map = 0;']);
        end
        d = 0;
        STStamp = datenum(Info.StartDate);
        ETStamp = datenum(Info.EndDate);
        for TS = STStamp:ETStamp  
            d = d+1;
            dirname = [Info.ModelStatesDir filesep 'LSM States' filesep 'Daily' filesep Info.NameIdentifier];
            for v = 1:numel(MVars)
                if MVars(v).map
                    disp(['Getting ' MVars(v).Name ' data for ' datestr(TS,'yyyy-mm-dd') ' - ' Info.NameIdentifier]);
                    eval(['load([dirname filesep datestr(TS,''yyyymmdd'')],''' MVars(v).Name '_output'');']);
                    eval([MVars(v).Name '_day(isnan(' MVars(v).Name '_output)) = 0;']);
                    if strcmp(Info.outputType,'avg') || strcmp(Info.outputType,'sum')
                        eval([MVars(v).Name '_map = ' MVars(v).Name '_map + ' MVars(v).Name '_output;']);
                    elseif strcmp(Info.outputType,'max')
                        eval([MVars(v).Name '_map = max(' MVars(v).Name '_map, ' MVars(v).Name '_output);']);
                    end
                end
            end
        end

        for v = 1:numel(MVars)
            if MVars(v).map 
                if strcmp(Info.outputType,'avg')
                    eval([MVars(v).Name '_map = ' MVars(v).Name '_map / d;']);
                end

                eval(['val = ' MVars(v).Name '_map;']);  	
                fname = [Info.outputType '_' datestr(STStamp,'yyyymmdd') '_' datestr(ETStamp,'yyyymmdd')];
                eval(['fpath = [Info.DisplayDir filesep ''GIS'' filesep Info.NameIdentifier filesep ''' MVars(v).Name '_composite''];']);
                print_maps(Info,val,Z+veght,fpath,fname);
            end
        end
    end
end

            
function [rain_t,snow_t,airt_t,pres_t,rh_t,tdew_t,vapp_t,ugrd_t,vgrd_t,dswrf_t,dlwrf_t] = get_daily_forcing_data(gen,Info,FVars,modelpars,SpatialData,Z,veght,TS)
% Gets the forcing data, or outputs daily data, and optionally makes hourly maps

% Initialize all values (these are for daily sums and averages)
rain_t = 0;
snow_t = 0;
airt_t = 0;
pres_t = 0;
rh_t = 0;
tdew_t = 0;
vapp_t = 0;
wind_t = 0;
ugrd_t = 0;
vgrd_t = 0;
dswrf_t = 0;
dlwrf_t = 0;

% Get hourly forcing variables
for j = 1:gen.DAY/gen.TS
    TS_sub = TS+((j-1)/(gen.DAY/gen.TS));
    % Only get the ones that are requested
    ForcingVars = [0 0 0 0 0 0 0 0];
   	if FVars(1).map == 1, ForcingVars(6) = 1; ForcingVars(2) = 1; end
    if FVars(2).map == 1, ForcingVars(6) = 1; ForcingVars(2) = 1; end
    if FVars(3).map == 1, ForcingVars(2) = 1; end
    if FVars(4).map == 1, ForcingVars(1) = 1; end
    if FVars(5).map == 1, ForcingVars(1) = 1; ForcingVars(2) = 1; ForcingVars(5) = 1; end
    if FVars(6).map == 1, ForcingVars(1) = 1; ForcingVars(2) = 1; ForcingVars(5) = 1; end
    if FVars(7).map == 1, ForcingVars(1) = 1; ForcingVars(2) = 1; ForcingVars(5) = 1; end
    if FVars(8).map == 1, ForcingVars(3) = 1; ForcingVars(4) = 1; end
    if FVars(9).map == 1, ForcingVars(3) = 1; ForcingVars(4) = 1; end
    if FVars(10).map == 1, ForcingVars(8) = 1; end
    if FVars(11).map == 1, ForcingVars(7) = 1; end
    % Get forcing data for each timestep
    [prec,airt,pres,vapp,ugrd,vgrd,dswrf,dlwrf] = get_forcing_data(Info,SpatialData,TS_sub,ForcingVars);
    % Right now, all adjustment factors are 1 (for multiplication) and 0
    % (for addition)
    rain_fraction = max(0,min(1,(1 + airt - (modelpars.tsnow + gen.KELVIN))/2));
    esat = 0.6108*exp(17.27*(airt-gen.KELVIN)./(237.3+(airt-gen.KELVIN)))*1000; % in Pa
    rh = vapp ./ esat;
    tdew = (log(vapp/1000)+0.49299)./(0.0707-0.00421*log(vapp/1000)); % [deg-C]
    rain = max(0,prec .* rain_fraction .* modelpars.RainMult);
    snow = max(0,prec .* (1-rain_fraction) .* modelpars.SnowMult);
    airt = airt + modelpars.TempAdjustment;
    dswrf = dswrf .* modelpars.SRadMult;
    dlwrf = dlwrf .* modelpars.LRadMult;
    
    % If specified, output hourly maps
    for v = 1:numel(FVars)
        if FVars(v).map && strcmp(Info.outputType,'ts')
            eval(['val = ' lower(FVars(v).Name) ';']);  
            disp(['Getting ' FVars(v).Name ' data for ' datestr(TS_sub,'yyyy-mm-dd HH:MM') ' - ' Info.NameIdentifier]);
            eval(['fname = datestr(TS_sub,''yyyymmddHHMM'');']);
            eval(['fpath = [Info.DisplayDir filesep ''GIS'' filesep Info.NameIdentifier filesep ''' FVars(v).Name '_ts''];']);
            print_maps(Info,val,Z+veght,fpath,fname);
        end
    end

    % Sum over the day ...
    rain_t = rain_t + rain;
    snow_t = snow_t + snow;
    airt_t = airt_t + airt;
    pres_t = pres_t + pres;
    rh_t = rh_t + rh;
    tdew_t = tdew_t + tdew;
    vapp_t = vapp_t + vapp;
    ugrd_t = ugrd_t + ugrd;
    vgrd_t = vgrd_t + vgrd;
    dswrf_t = dswrf_t + dswrf;
    dlwrf_t = dlwrf_t + dlwrf;
end
% and divide by 24 hours to get the average
airt_t = airt_t / (gen.DAY/gen.TS);
pres_t = pres_t / (gen.DAY/gen.TS);
rh_t = rh_t / (gen.DAY/gen.TS);
tdew_t = tdew_t / (gen.DAY/gen.TS);
vapp_t = vapp_t / (gen.DAY/gen.TS);
ugrd_t = ugrd_t / (gen.DAY/gen.TS);
vgrd_t = vgrd_t / (gen.DAY/gen.TS);
dswrf_t = dswrf_t / (gen.DAY/gen.TS);
dlwrf_t = dlwrf_t / (gen.DAY/gen.TS);

function [r] = getsolarindexes(Info,SpatialData,modelpars,includeveg)
    % Function to get the solar forcing indexes

    STStamp = datenum(Info.StartDate);
    ETStamp = datenum(Info.EndDate);
    veght = SpatialData.canopy_height;
    % For now, vegcover_canopy is constant
    vegcover_canopy = SpatialData.canopy_coverage/100;
    gen.TINY = 1E-5;

    % For Vegetation Dynamics model
    dh_LAI_c = 1;
    dx_LAI_c = modelpars.LAI_c_tmax-modelpars.LAI_c_tmin;
    % V0.6: We don't use lai of ground vegetation here
    r_0 = 0;
    r = 0;

    for TS = STStamp:ETStamp
        disp(['Computing solar forcing index for ' Info.NameIdentifier ' - ' datestr(TS)]);
        Doy = date2doy(TS);
        % Get solar forcing for flat ground
        [R_0,Id,~] = solarradiation_direct(0,0,SpatialData.LatMap,Doy,1/2,24-1/2,1,single([]),uint8([]));
        Ir = 1-Id;
        Diffuse_fr = 1 - Ir;
        R_0_diffuse = R_0.*Diffuse_fr;
        r_0 = r_0 + R_0;

        % Fo now, no vegetation dynamics
        % TODO: Add Dynamics based on temperature here
        state.AvMaxT = 15;
        Alpha_c = max(0,min(1,(dh_LAI_c/dx_LAI_c)*(state.AvMaxT-modelpars.LAI_c_tmin) - dh_LAI_c/(2*pi)*sin((2*pi/dx_LAI_c)*(state.AvMaxT-modelpars.LAI_c_tmin))));
        LAI_0_c = (Alpha_c * modelpars.maxLAI_c + (1-Alpha_c) * modelpars.minLAI_c);

        if includeveg == 0
            % In case of bare earth solar forcing index, compute solar forcing
            % with topography alone
            R = solarradiation_direct(atan(SpatialData.S),SpatialData.R,SpatialData.LatMap,Doy,1/2,24-1/2,1,SpatialData.Alphas,SpatialData.AlphaMat);
            r = r + R;
        elseif includeveg == 1
            % In case of veg solar forcing index, compute solar forcing with
            % veg + topography

            % Compute the skyview factor
            SkyView_canopy = 0;
            SkyView_noveg = 0;
            for i = 1:numel(SpatialData.Alphas)-1
                SkyView_noveg = SkyView_noveg + (pi/2 - (double(SpatialData.AlphaMat(:,:,i))/(255/(pi/2)))) / (pi/2) / (numel(SpatialData.Alphas)-1); 
                SkyView_canopy = SkyView_canopy + (pi/2 - (double(SpatialData.VegAlphaMat(:,:,i))/(255/(pi/2)))) / (pi/2) / (numel(SpatialData.Alphas)-1); 
            end
            adj_cover = max(max(0,min(1,veght/15)),vegcover_canopy);
            SkyView = min(SkyView_canopy,SkyView_noveg) .* (1-adj_cover).^2;

            % Compute tau for diffuse shortwave radiation
            G = 0.5;
            Alpha = 0.5;
            K_p = (1-Alpha).^1/2;
            LAI_c = vegcover_canopy * LAI_0_c;
            KGL = K_p .* G .* LAI_c .* (1-SkyView);
            tau_b_shortwave = ((1-KGL) .* exp(-KGL)) + ((KGL).^2 .* min(1,expint(KGL)));
            tau_b_shortwave(isnan(tau_b_shortwave)) = 0; tau_b_shortwave(tau_b_shortwave<0) = 0; tau_b_shortwave(tau_b_shortwave>1) = 1;    

            Rs_d_direct = 0;
            Rs_d_diffuse = 0;
            for hour = 1:24
                % Get bare earth solar forcing and solar forcing assuming veg
                % elements are opaque, we will compare these two to find out
                % how much vegetation is between the sun and each pixel (this
                % gives vegetation cover in direction of sun)
                [R_NoVeg,~,alpha] = solarradiation_direct(atan(SpatialData.S),SpatialData.R,SpatialData.LatMap,Doy,hour-1/2,hour+1/2,1,SpatialData.Alphas,SpatialData.AlphaMat);
                R_Top = solarradiation_direct(atan(SpatialData.S_veg),SpatialData.R_veg,SpatialData.LatMap,Doy,hour-1/2,hour+1/2,1,SpatialData.Alphas,SpatialData.VegAlphaMat);
                R_Shadows = R_Top;      % V0.6: Remove Dependence on VegCover
                % Get direct solar radiation only
                Ir = 1-Id;
                Diffuse_fr = 1 - Ir;
                R_direct = R_NoVeg.*(1-Diffuse_fr);
                R_diffuse = R_NoVeg.*Diffuse_fr;

                Alpha = 0.5;                                % Transmittance
                G = 0.5;                                	% Geometry factor
                K_p = (1-Alpha).^1/2;
                alpha(alpha<gen.TINY) = gen.TINY;
                K_b = G./alpha;
                % This is vegetation cover in direction of sun
                adj_cover = max(max(0,min(1,veght/15)),vegcover_canopy);           % V0.6: Fix occasional visual artifacts around the edges of trees caused by mismatch of lidar cover and height
                vegcover_sun = (1-min(1,R_Shadows./R_NoVeg).*(1-adj_cover).^2);    % Vegetation cover between the sun and each pixel
                % Effective leaf area index in direction of sun
                LAI_cs = LAI_0_c .* vegcover_sun;
                LAI_cs(isnan(LAI_cs)) = 0;
                % Tau for direct radiation
                tau_a = exp(-K_p .* K_b .* LAI_cs); 
                tau_a(isnan(tau_a)) = 0; tau_a(tau_a<0) = 0; tau_a(tau_a>1) = 1;

                Rs_d_direct = Rs_d_direct + (R_direct .* tau_a); 
                Rs_d_diffuse = Rs_d_diffuse + (R_diffuse .* tau_b_shortwave);

            end
            % Get daily average radiation
            Rs_d_direct = Rs_d_direct/24;
            Rs_d_diffuse = Rs_d_diffuse/24;
            r = r + Rs_d_direct+Rs_d_diffuse;  
        end
    end
    % Compare to flat surface solar forcing
    r = r ./ r_0;


function print_maps(Info,val,Z,fpath,fname)
    % Function to create GIS Files

    % Area outside of model domain -> nodata
    val(isnan(Z)) = NaN;

    if ~exist(fpath,'file')
        mkdir(fpath)
    end
    % Make sure that output folder exists, write to tiff file
    if strcmp(Info.outputFType,'tif')
        ofname = [fpath filesep fname '.tif'];
        if ~exist(ofname,'file')
            ifname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'DEM.tif'];
            imwrite2tif(val,[],ofname,'single');
            eval(['!python "' Info.ProgramFilesDir filesep 'Python' filesep 'gdalcopyproj.py" "' ifname '" "' ofname '"']);
        end
    elseif strcmp(Info.outputFType,'mat')
        ofname = [fpath filesep fname '.mat'];
        save(ofname,'-v6','val');
    end
