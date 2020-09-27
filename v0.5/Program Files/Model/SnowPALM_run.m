function SnowPALM_run(gen,Info,SpatialData,modelpars,MVars)   
%   Function to run the core of SnowPALM.  Included in the following
%   function is the model physics related to vegetation dynamics,
%   radiation, interception, snow interception, and snowpack accounting. 
%   The entire model is fully distributed and runs at hourly resolution. 
%
% Inputs: 'Info' is a structure that contains program parameters
%         'modelpars' contain model parameters (see 'get_model_pars.m' for more information)
%         'SpatialData' is a structure that contains all necissary spatial data
%         'gen.TS' is the model timestep [s]
%         'ForcingData' is a structure that contains all necissary forcing data
% No Outputs
%
%   Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated November 2012
%

%% Set Up Model
Info.ProgDispFlag       = 1;            % Display Program Progress
STStamp = datenum(Info.StartDate);      % Start Date
ETStamp = datenum(Info.EndDate);        % End Date 

% Set up the output data structures
if isfield(SpatialData,'DispLocs')
    Info.ExportCSV = 1; 
else
    Info.ExportCSV = 0; 
end

% Update Spatial data structure to reflect that we only do computations on 
% the model domain, not the whole tile if not complete coverage
% TODO: Not all of these inputs are necissary (figure out which ones are
% needed)
sz = size(SpatialData.z);
nanlocs = isnan(SpatialData.z);
SpatialData.z = SpatialData.z(~nanlocs);
SpatialData.x = SpatialData.x(~nanlocs);
SpatialData.y = SpatialData.y(~nanlocs);
SpatialData.S = SpatialData.S(~nanlocs);
SpatialData.R = SpatialData.R(~nanlocs);
SpatialData.LonMap = SpatialData.LonMap(~nanlocs);
SpatialData.LatMap = SpatialData.LatMap(~nanlocs);
SpatialData.canopy_height = SpatialData.canopy_height(~nanlocs);
SpatialData.canopy_coverage = SpatialData.canopy_coverage(~nanlocs);
for i = 1:numel(SpatialData.Alphas)-1
    tmp = SpatialData.AlphaMat(:,:,i);
    AlphaMat(:,i) = tmp(~nanlocs);
end
SpatialData.AlphaMat = AlphaMat;
SpatialData.CanPos = SpatialData.CanPos(~nanlocs);
SpatialData.WgtFun = SpatialData.WgtFun(~nanlocs);
SpatialData.WgtFun_shortwave = SpatialData.WgtFun_shortwave(~nanlocs);
SpatialData.S_0 = SpatialData.S_0(~nanlocs);
SpatialData.R_0 = SpatialData.R_0(~nanlocs);
SpatialData.S_veg = SpatialData.S_veg(~nanlocs);
SpatialData.R_veg = SpatialData.R_veg(~nanlocs);
for i = 1:numel(SpatialData.VegAlphas)-1
    tmp = SpatialData.VegAlphaMat(:,:,i);
    VegAlphaMat(:,i) = tmp(~nanlocs);
end
SpatialData.VegAlphaMat = VegAlphaMat;
SpatialData.Sx = SpatialData.Sx(~nanlocs);
SpatialData.Sb = SpatialData.Sb(~nanlocs);
if Info.ExportCSV
    for i = 1:numel(SpatialData.DispLocs)
        pos = zeros(sz);
        weights = zeros(sz);
        pos(SpatialData.DispLocs(i).pos) = 1;
        weights(SpatialData.DispLocs(i).pos) = SpatialData.DispLocs(i).weights;
        pos = pos(~nanlocs);
        weights = weights(~nanlocs);
        SpatialData.DispLocs(i).pos = find(pos);
        SpatialData.DispLocs(i).weights = weights(pos == 1);
    end
end
sz_domain = size(SpatialData.z);

% If using distributed forcing correction factors
if SpatialData.POIForcingCorrectionFactors
    TempAdjustment = 0;
    TSnow = 0;
    PrecipMult = 0;
    TotalWeights = 0;
    for i = 1:numel(SpatialData.DispLocs)
        if ~isempty(SpatialData.DispLocs(i).PrecipMult)
            dist = sqrt((SpatialData.x-SpatialData.DispLocs(i).X).^2 + (SpatialData.y-SpatialData.DispLocs(i).Y).^2);
            weight = 1./(dist.^Info.forcingCorrExponent); 
            PrecipMult = PrecipMult + weight * SpatialData.DispLocs(i).PrecipMult;
            TempAdjustment = TempAdjustment + weight * SpatialData.DispLocs(i).TempAdjustment;
            TSnow = TSnow + weight * SpatialData.DispLocs(i).TSnow;
            TotalWeights = TotalWeights + weight;
        end
    end
    modelpars.PrecipMult = PrecipMult ./ TotalWeights;
    modelpars.TempAdjustment = TempAdjustment ./ TotalWeights;
    modelpars.tsnow = TSnow ./ TotalWeights;
    modelpars.RainMult = modelpars.PrecipMult;
    modelpars.SnowMult = modelpars.PrecipMult;
end

% Set up structure needed for the model forcing data 
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
    ForcingVars = [1 1 1 1 1 1 1 1];
end
% Used when constructing file name to read from
for j = 1:numel(Info.Forcing.var)
    ForcingDir = Info.ForcingDirs{Info.Forcing.source{j}};
    fprefix = strsplit(ForcingDir,filesep);
    fprefix = strsplit(char(fprefix(end)),'/');
    Info.Forcing.fprefix{j} = char(fprefix(end));
end

% Name of model restart files
restartdirname = [Info.ModelStatesDir filesep 'LSM States' filesep 'Restart' filesep Info.NameIdentifier]; 

% For Tabular Data
TabularData = [];
if Info.ExportCSV 
    cs = SpatialData.cellsize;
    for i = 1:numel(SpatialData.DispLocs)
        TabularData(i).Name = SpatialData.DispLocs(i).Name;
        TabularData(i).area = numel(SpatialData.DispLocs(i).pos(:)) * cs ^ 2;
      	TabularData(i).TS = STStamp:gen.TS/gen.DAY:ETStamp+1-(gen.TS/gen.DAY);
        for v = 1:numel(MVars)
            if ~isempty(SpatialData.DispLocs(i).pos)
                eval(['TabularData(i).' MVars(v).Name ' = NaN(numel(STStamp:ETStamp)*gen.DAY/gen.TS,1);']);
            end
        end
    end
end

% For Spatial Data
for v = 1:numel(MVars)
    if MVars(v).map
        if strcmp(Info.MapType,'day')
            eval([MVars(v).Name  '_output = single(zeros(sz));']);
            outdirname = [Info.ModelStatesDir filesep 'LSM States' filesep 'Daily' filesep Info.NameIdentifier];
            if ~exist(outdirname,'file')
                mkdir(outdirname); 
            end
        elseif strcmp(Info.MapType,'ts')
            eval([MVars(v).Name  '_output = single(zeros([sz gen.DAY/gen.TS]));']);
            outdirname = [Info.ModelStatesDir filesep 'LSM States' filesep 'TS' filesep Info.NameIdentifier];
            if ~exist(outdirname,'file')
                mkdir(outdirname); 
            end
        end
    end
end

% Initialize model states (initialize snow states to zero)
state.Q_snow = zeros(sz_domain);
state.swe = zeros(sz_domain);
state.liq_swe = zeros(sz_domain);
state.sdepth = zeros(sz_domain);
state.T_snows = zeros(sz_domain);
state.T_snowm = zeros(sz_domain);
state.T_snowb = zeros(sz_domain);
state.T_soil = zeros(sz_domain);
state.Q_soil = zeros(sz_domain);
state.ice_percent_soil = zeros(sz_domain);
state.cansnowstor = zeros(sz_domain);
state.canstor = zeros(sz_domain);
state.dcorr = zeros(sz_domain);
state.acorr = zeros(sz_domain);
state.swe_age_a = zeros(sz_domain);
state.swe_age_d = zeros(sz_domain);
state.smc = ones(sz_domain) .* modelpars.initialSMC;
state.Temp_c = [];
state.AvMaxT = [];

% Prevent zeros
modelpars.maxLAI_c = max(gen.TINY,modelpars.maxLAI_c);
modelpars.minLAI_c = max(gen.TINY,modelpars.minLAI_c);

% Convert canopy cover data (0..100 -> 0..1)
cover = SpatialData.canopy_coverage/100;
veght = SpatialData.canopy_height;
    
% For Vegetation Dynamics model
dh_LAI_c = 1;
dx_LAI_c = modelpars.LAI_c_tmax-modelpars.LAI_c_tmin;
dh_LAI_g = 1;
dx_LAI_g = modelpars.LAI_g_tmax-modelpars.LAI_g_tmin;

% For wind parameterization
if Info.UseWindModel
    sx = SpatialData.Sx;
    sb = SpatialData.Sb * modelpars.Sb_Multiplier + 1;
end

if Info.InitializeFromRestartFile
    TS = STStamp;
    restartfilename = [restartdirname filesep datestr(TS,'yyyymmdd') '.mat'];
    if exist(restartfilename,'file')
        load(restartfilename)
    end
end     

%% Main model Time Loop
tic
pause(1.1)
f = 0;  % Timesteps from the beginning of the model run
for TS = STStamp:ETStamp
    
    % Get model state from previous timestep if specified
    restartfilename = [restartdirname filesep datestr(TS,'yyyymmdd') '.mat'];
    restartfilename_next = [restartdirname filesep datestr(TS+1,'yyyymmdd') '.mat'];
    modeloutputfile = [outdirname filesep datestr(TS,'yyyymmdd') '.mat'];
    
    if Info.ContinuationRun
        if exist(modeloutputfile,'file')
            exec_model = 0;
        elseif exist(restartfilename,'file')
            load(restartfilename)
            TabularData = TabularData_disp;
            exec_model = 1;
        else
            warning(['Cannot Complete Continuation Run (' restartfilename ' does not exist!']);
            exec_model = 1;
        end
    else
        exec_model = 1;
    end
    
    if exec_model
        TS_string = datestr(TS,'yyyymmdd');
        Doy = date2doy(TS); 

        % If the program progress is to be displayed, print out every second
        if Info.ProgDispFlag
            tmp = toc;
            if tmp > 1
                disp(['Running Land-Surface Model for ' Info.NameIdentifier ': ' TS_string]);
                tic
            end
        end

        % Keep track of the daily maximum and minimum daily air temperatures
        minairt = ones(sz_domain)*500;
        maxairt = ones(sz_domain);

        %% Get Forcing Data for the day
        [ForcingVals] = get_daily_nldas_forcing(Info,SpatialData,TS,ForcingVars);
        if isempty(state.AvMaxT)
            state.AvMaxT = ForcingVals(1).airt + modelpars.TempAdjustment;
        end
        if isempty(state.Temp_c)
            state.Temp_c = ForcingVals(1).airt + modelpars.TempAdjustment;
        end

        %% Vegetation Dynamics

        % Vegetation height and cover of both canopy and understory (cover is
        % scaled with LAI (next section)
        vegcover_understory_0 = 1-max(0,min(1,((state.sdepth/gen.M2MM)./modelpars.groundveght)));
        vegcover_canopy_0 = cover;
        veght_understory = max(0,modelpars.groundveght - (state.sdepth/gen.M2MM));
        veght_canopy = max(0,veght - (state.sdepth/gen.M2MM));
        vegcover_canopy_0(veght_canopy==0) = 0;
        vegcover_tot_0 = min(1,vegcover_understory_0+vegcover_canopy_0);
        veght_tot = max(veght_canopy,veght_understory);

        % Compute seasonally varying LAI
        Alpha_c = max(0,min(1,(dh_LAI_c/dx_LAI_c)*(state.AvMaxT-modelpars.LAI_c_tmin) - dh_LAI_c/(2*pi)*sin((2*pi/dx_LAI_c)*(state.AvMaxT-modelpars.LAI_c_tmin))));
        Alpha_g = max(0,min(1,(dh_LAI_g/dx_LAI_g)*(state.AvMaxT-modelpars.LAI_g_tmin) - dh_LAI_g/(2*pi)*sin((2*pi/dx_LAI_g)*(state.AvMaxT-modelpars.LAI_g_tmin))));
        LAI_0_c = (Alpha_c * modelpars.maxLAI_c + (1-Alpha_c) * modelpars.minLAI_c);
        LAI_c = vegcover_canopy_0 .* LAI_0_c;
        LAI_0_g = (Alpha_g * modelpars.maxLAI_g + (1-Alpha_g) * modelpars.minLAI_g);
        LAI_g = vegcover_understory_0 .* LAI_0_g;
        vegcover_understory = vegcover_understory_0 .* (LAI_0_c ./ modelpars.maxLAI_c);
        vegcover_canopy = vegcover_canopy_0 .* (LAI_0_c ./ modelpars.maxLAI_c);

        windlevel = modelpars.windlevel + veght_canopy;
        ka_par_s = modelpars.karman^2 ./ ( log( ( windlevel ) / modelpars.sroughness) ).^2;
        cansnowstorcap = 4.4*LAI_c;     % Canopy snow storage

        %% Skyview Factor
        SkyView_canopy = 0;
        SkyView_noveg = 0;
        for i = 1:numel(SpatialData.Alphas)-1
            SkyView_noveg = SkyView_noveg + (pi/2 - (double(SpatialData.AlphaMat(:,i))/(255/(pi/2)))) / (pi/2) / (numel(SpatialData.Alphas)-1); 
            SkyView_canopy = SkyView_canopy + (pi/2 - (double(SpatialData.VegAlphaMat(:,i))/(255/(pi/2)))) / (pi/2) / (numel(SpatialData.Alphas)-1); 
        end
        adj_cover = max(max(0,min(1,veght/15)),vegcover_canopy);
        SkyView = min(SkyView_canopy,SkyView_noveg) .* (1-adj_cover).^2;

        
    %     var = NaN(sz);
    %     var(~nanlocs) = SkyView;
    %     pcolor(flipud(double(var)))
    %     shading flat
    %     colorbar
    %     pause

        %% Tau for longwave and diffuse shortwave radiation

        G = 0.5;
        Alpha = 0.5;
        K_p = (1-Alpha).^1/2;
        KGL = K_p .* G .* LAI_0_c .* (1-SkyView);
        tau_b_shortwave = ((1-KGL) .* exp(-KGL)) + ((KGL).^2 .* min(1,expint(KGL)));
        tau_b_shortwave(isnan(tau_b_shortwave)) = 0; tau_b_shortwave(tau_b_shortwave<0) = 0; tau_b_shortwave(tau_b_shortwave>1) = 1;

        G = 0.5;  
        Alpha = 0.5;
        K_p = (1-Alpha).^1/2;
        KGL = K_p .* G .* (1-SkyView) .* LAI_0_c;
        tau_b_longwave = ((1-KGL) .* exp(-KGL)) + ((KGL).^2 .* min(1,expint(KGL)));
        tau_b_longwave(isnan(tau_b_longwave)) = 0; tau_b_longwave(tau_b_longwave<0) = 0; tau_b_longwave(tau_b_longwave>1) = 1;

        %% Sub daily timestep
        for j = 1:gen.DAY/gen.TS;  
            f = f+1;

            %% Unpack Forcing Data
            airt = ForcingVals(j).airt + modelpars.TempAdjustment;
            vapp = ForcingVals(j).vapp;
            prec = ForcingVals(j).prec;
            rain_fraction = max(0,min(1,(1 + airt - (modelpars.tsnow + gen.KELVIN))/2));
            esat = 0.6108*exp(17.27*(airt-gen.KELVIN)./(237.3+(airt-gen.KELVIN)))*1000; % in Pa
            rh = vapp ./ esat;
            tdew = (log(vapp/1000)+0.49299)./(0.0707-0.00421*log(vapp/1000)); % [deg-C]
            rain = max(0,prec .* rain_fraction .* modelpars.RainMult);
            snow_0 = max(0,prec .* (1-rain_fraction) .* modelpars.SnowMult);
            pres = ForcingVals(j).pres;
            ugrd = ForcingVals(j).ugrd;
            vgrd = ForcingVals(j).vgrd;
            dswrf = ForcingVals(j).dswrf .* modelpars.SRadMult;
            dlwrf = ForcingVals(j).dlwrf .* modelpars.LRadMult;
            wind = sqrt(ugrd.^2+vgrd.^2);
            windangle = atan2(vgrd,ugrd);
            windangle = windangle + pi;
            windangle_p = windangle + 2*pi;
            windangle_p(windangle > 0) = 0;
            windangle(windangle < 0) = 0;
            windangle = windangle + windangle_p;
            minairt = min(minairt,airt);
            maxairt = max(maxairt,airt);
            % Adjust snowfall for wind influences
            if Info.UseWindModel
                snow = snow_0 .* sx .* sb;
            else
                snow = snow_0;
            end

            liq_swe_cap = state.swe .* modelpars.liq_swe_cap_mult;
            snowcover = (state.swe>0);

            %% Snow Depth
            new_depth = max(gen.TINY,snow/modelpars.density_i);
            new_frac = new_depth./(state.sdepth + new_depth);
            state.swe_age_d = min(state.swe_age_d + 1,modelpars.maxage_d*(gen.DAY/gen.TS)); 
            state.swe_age_d = state.swe_age_d - (new_frac .* state.swe_age_d);
            snowdensity = modelpars.density_i + (modelpars.maxdensity - modelpars.density_i) * state.swe_age_d/(modelpars.maxage_d*(gen.DAY/gen.TS));
            state.dcorr = (1-new_frac) .* state.dcorr;
            ripe_d_par = state.liq_swe ./ liq_swe_cap .* modelpars.ripe_d_par;
            state.dcorr = min(modelpars.maxdensity-snowdensity,state.dcorr + ripe_d_par * (gen.TS/gen.DAY)); 
            state.dcorr(state.swe == 0) = 0;
            snowdensity = snowdensity + state.dcorr;
            snowdensity_corr = min(modelpars.maxdensity,snowdensity);
            state.sdepth = (state.swe + state.liq_swe) ./ snowdensity_corr;   

            %% Albedos, Emmisivities, Roughness

            % Snowpack Albedo
            state.swe_age_a = state.swe_age_a + (gen.TS/gen.DAY); 
            state.swe_age_a = state.swe_age_a .* max(0,(modelpars.albedo_snow_reset-snow)/modelpars.albedo_snow_reset);
            state.acorr = state.acorr .* max(0,(modelpars.albedo_snow_reset-snow)/modelpars.albedo_snow_reset);
            albedosnow = (modelpars.albedo_i-modelpars.minalbedo).*exp(-1./(modelpars.maxage_a) .*state.swe_age_a) + modelpars.minalbedo;
            ripe_a_par = state.liq_swe ./ liq_swe_cap .* modelpars.ripe_a_par;
            state.acorr = min(albedosnow-modelpars.minalbedo,state.acorr + ripe_a_par * (gen.TS/gen.DAY)); 
            state.acorr(state.swe == 0) = 0;
            albedosnow = albedosnow - state.acorr;

            % Ground Albedo
            albedo_g = albedosnow.*snowcover + modelpars.albedo_ground.*(1-snowcover);

            % Canopy Albedo
            pcsnow = state.cansnowstor ./ cansnowstorcap;   
            pcsnow(isnan(pcsnow)) = 0;
            pcsnow(isinf(pcsnow)) = 0;
            albedo_c = (1-pcsnow) .* modelpars.albedo_tree + pcsnow .* modelpars.albedo_i; 

            % Emmissivities of the ground and canopy
            emiss_c = modelpars.emiss_tree;
            emiss_g = modelpars.emiss_snow.*snowcover + modelpars.emiss_ground.*(1-snowcover); 

            % Zero Plane displacement height and roughness length
            vzeroplan = 0.067 *veght_tot;
            vroughness = max(modelpars.sroughness,0.123 *veght_tot);

            %% Radiation

            if max(dswrf(:)) > 5
                tm = j*gen.TS/gen.HOUR + mean(SpatialData.LonMap(:))/180*12;
                if tm < 0
                    tm = 24 + tm; 
                end

                % Get bare earth solar forcing and solar forcing assuming veg
                % elements are opaque, we will compare these two to find out
                % how much vegetation is between the sun and each pixel (this
                % gives vegetation cover in direction of sun)
                [R_0,Id,~] = solarradiation_direct(0,0,SpatialData.LatMap,Doy,tm,tm,1,single([]),uint8([]));
                [R_NoVeg,~,alpha] = solarradiation_direct(atan(SpatialData.S),SpatialData.R,SpatialData.LatMap,Doy,tm,tm,1,SpatialData.Alphas,SpatialData.AlphaMat);
                R_Top = solarradiation_direct(atan(SpatialData.S_veg),SpatialData.R_veg,SpatialData.LatMap,Doy,tm,tm,1,SpatialData.Alphas,SpatialData.VegAlphaMat);
                R_Shadows = R_Top .* (1-vegcover_canopy);
                % Get direct solar radiation only
                Ir = 1-Id;
                Diffuse_fr = 1 - Ir;

                R_NoVeg = R_NoVeg + R_0.*Diffuse_fr;
                R_Top = R_Top + R_0.*Diffuse_fr;
                R_Shadows = R_Shadows + R_0.*Diffuse_fr;
                SFI_Top = R_Top./R_0; 
                SFI_Top = min(SFI_Top,5);
                SFI_NoVeg = R_NoVeg ./ R_0;
                SFI_NoVeg = min(SFI_NoVeg,5);

                Rs_d_direct_top = max(0,dswrf .* SFI_Top .* (1-Diffuse_fr));
                Rs_d_direct = max(0,dswrf .* SFI_NoVeg .* (1-Diffuse_fr));
                Rs_d_diffuse = max(0,dswrf .* Diffuse_fr);

                Alpha = 0.5;                                % Transmittance
                G = 0.5;                                	% Geometry factor
                K_p = (1-Alpha).^1/2;
                alpha(alpha<gen.TINY) = gen.TINY;
                K_b = G./alpha;
                % This is vegetation cover in direction of sun
                vegcover_sun = (1-min(1,R_Shadows./R_NoVeg).*(1-vegcover_canopy));    % Vegetation cover between the sun and each pixel
                % Effective leaf area index in direction of sun
                LAI_cs = LAI_0_c .* vegcover_sun;
                LAI_cs(isnan(LAI_cs)) = 0;
                % Tau for direct radiation
                tau_a = exp(-K_p .* K_b .* LAI_cs); 
                tau_a(isnan(tau_a)) = 0; tau_a(tau_a<0) = 0; tau_a(tau_a>1) = 1;

                % Transmitted direct radiation through the canopy

                Rs_g_direct = Rs_d_direct .* tau_a;
                Rs_g_diffuse = Rs_d_diffuse .* tau_b_shortwave;
                Rs_down_g = max(0,Rs_g_direct + Rs_g_diffuse);
                Rs_down_c = (Rs_d_direct_top + Rs_d_diffuse - Rs_down_g);

    %             var = NaN(sz);
    %             var(~nanlocs) = Rs_g_direct;
    %             pcolor(flipud(double(var)))
    %             shading flat
    %             colorbar
    %             pause

                albedo_g_r = min(1,albedo_g+max(0,(1-alpha) * 0.5));
                albedo_c = min(1,albedo_c+max(0,(1-alpha) * 0.5));

                Rs_scatter = Rs_down_g .* albedo_g_r .* (modelpars.ScatterEff);
                R_s_up_g = Rs_down_g .* albedo_g_r .* (1-modelpars.ScatterEff);
                Rs_g = Rs_down_g - R_s_up_g;
                R_s_up_c = Rs_down_c .* albedo_c;
                Rs_c = Rs_down_c - R_s_up_c + R_s_up_g .* albedo_c .* (1-tau_a);

                Edgeness_Par = max(0,min(1,SFI_Top .^ (exp(-veght./modelpars.EFLength_vert))-1));
                ShadowIndex = vegcover_sun.* exp(-veght./modelpars.EFLength_vert);
            else 
                ShadowIndex = zeros(sz_domain); 
                Rs_scatter = zeros(sz_domain); 
                Rs_d_direct_top = zeros(sz_domain);
                Rs_d_diffuse = zeros(sz_domain);
                Rs_g = zeros(sz_domain);
                Rs_c = zeros(sz_domain);
                Edgeness_Par = zeros(sz_domain);
            end

    %      	%% Canopy Liquid Storage (off for now)
    %         
    %         % Account for canopy storage of rain and snow
    %         % Treat canopy storages like buckets that fill to a capacity, and
    %         % either drip to the ground or evaporate or sublimate
    %         
    %         canstorcap = 0.2 * LAI_0_c;  % CanStorcap [mm]   % Dickinson, 1984
    % 
    %         tfall = (1-cover) .* rain;      	% find rain that falls on canopy
    %         % Now talk in terms of per unit area of canopy
    %         can_rain = rain;
    %         state.canstor = state.canstor + can_rain;                       % If canopy rain completely fills canopy, then it becomes throughfall
    %         tfall_canopy = max(0,state.canstor - canstorcap);
    %         state.canstor = min(canstorcap, state.canstor);
    %         
    %         drip = min(modelpars.candriprate * gen.TS/gen.HOUR,state.canstor); 	% calculate canopy drip rate
    %         state.canstor = max(0,state.canstor - drip);
    % 
    %         omega_wc = (state.canstor ./ canstorcap).^(2/3);
    %         omega_wc = max(omega_wc,0);
    %         omega_wc(pevap_c<0) = 1;                              % completely wet when dew occurs
    %         acevap = omega_wc .* pevap_c ;                        % actual wet canopy evaporation in mm    
    %         acevap(acevap<0)=0;
    %         state.canstor = max(0,state.canstor - acevap);                  % Subtract evaporated water from the canopy
    %         
    %         tfall_canopy = tfall_canopy + drip;
    %         tfall = tfall + vegcover_canopy .* tfall_canopy;                       % Account for fractional coverage of canopy
            tfall = rain;

            %% Canopy Snow Storage

            snow_c = snow .* cover;

            locs = find((state.cansnowstor > 0 | snow > 0) & LAI_c > 0);

            tsfall_canopy = zeros(sz_domain);
            melt_drip = zeros(sz_domain);
            snow_unload = zeros(sz_domain);
            acsub = zeros(sz_domain);

            tsfall_open = (1-cover) .* snow;                              % find snow that falls on canopy

            if ~isempty(locs)
                cstor_prev = state.cansnowstor;
                [tsfall_canopy(locs),melt_drip(locs),snow_unload(locs),state.cansnowstor(locs),acsub(locs)] = ...
                    getcsnowstor(gen,modelpars,state.Temp_c(locs),Rs_c(locs),LAI_c(locs),snow_c(locs),airt(locs),...
                    vapp(locs),state.cansnowstor(locs),wind(locs),rh(locs));

                acsub(isnan(acsub)) = 0;
                balance = max(0,cstor_prev+snow_c-tsfall_canopy-melt_drip-snow_unload-state.cansnowstor-acsub);
                tsfall_canopy = max(0,tsfall_canopy + balance);            % If snow in canopy gets burried by deepening pack
            end

            tsfall = tsfall_open + tsfall_canopy;

            % Keep track of rainfall, meteoric rain on snow, and rain+drip on snow
            rain_on_snow = rain;
            rain_on_snow(state.swe == 0) = 0;
            tfall = tfall + melt_drip;
            rain_on_snow_sfc = tfall;
            rain_on_snow_sfc(state.swe == 0) = 0;
            state.swe = max(0,state.swe + tsfall + snow_unload);  

            %% Skin Temperatures

            % Canopy Temperature
            ka = modelpars.karman^2 * wind ./ ( log( (modelpars.windlevel-vzeroplan) / modelpars.vroughness) ).^2;
            ka(isnan(ka)) = 0;

            C1 = (Rs_d_direct_top+Rs_d_diffuse).*(1-modelpars.albedo_tree) - (emiss_c .* modelpars.boltz .* airt.^4);
            C2 = emiss_c .* modelpars.boltz .* ((airt+5).^4-(airt-5).^4)/10;    % Linearized sigma-T curve
            C3 = modelpars.kappa_d_tree;
            Cx = 2;
            C4 = (modelpars.rhoa * modelpars.specheat_a .* ka + Cx);

            state.Temp_c = airt + C1./(C4+C3-C2);

            Temp_c_cansnowstor = min(gen.KELVIN,state.Temp_c);
            state.Temp_c(state.cansnowstor>0) = 0;
            Temp_c_cansnowstor(state.cansnowstor == 0) = 0;
            state.Temp_c = state.Temp_c + Temp_c_cansnowstor;

            Temp_c_factor = exp(-veght./modelpars.EFLength_vert);
            state.Temp_c = (state.Temp_c .* Temp_c_factor + airt .* (1-Temp_c_factor));
            sigmaT_c = emiss_c .* modelpars.boltz .* state.Temp_c.^4;       
            Rl_g = (dlwrf .* tau_b_longwave + sigmaT_c .* (1-tau_b_longwave));
            Rc_add_g = max(0,Edgeness_Par .* Rs_d_direct_top .* modelpars.EdgenessImportance);% + (Rs_g .* albedo_g .* modelpars.EdgenessImportance);
            Rnet_in_g = Rs_g + Rl_g + Rc_add_g.*ShadowIndex;

            T_soil = state.Q_soil./ ((state.smc .* modelpars.groundlayer_thickness * modelpars.specheat_i * modelpars.rhoi) + (1-state.smc .* modelpars.groundlayer_thickness) * modelpars.specheat_i * modelpars.rhoi);

            %% Energy Balance Snow Model

            locs = find(state.swe > 0 | tsfall > 0);
            melt = zeros(sz_domain);
            sublimation = zeros(sz_domain);
            H = zeros(sz_domain);
            lambdaE = zeros(sz_domain);
            ka = zeros(sz_domain);
            g = zeros(sz_domain);
            Q_m = zeros(sz_domain);
            Rn_snow = zeros(sz_domain);
            Q_soil_add = zeros(sz_domain);
            A_p = zeros(sz_domain);
            Q_conv = zeros(sz_domain);

            if ~isempty(locs)
                Ts0 = airt(locs);
                done = 0;
                niters = 20; thresh = 1;
                iter = 0;

                Ts = Ts0 - 2;
                dT = 4;
                [~,~,~,~,~,~,~,~,~,H0,lambdaE0,~,~,~,~,Rn_snow0] = get_snow(gen,Ts,T_soil(locs),state.swe(locs),state.liq_swe(locs),state.Q_snow(locs),...
                    state.Q_soil(locs),snowdensity(locs),Rnet_in_g(locs),rain_on_snow_sfc(locs),tsfall(locs),airt(locs),vapp(locs),wind(locs),pres(locs),modelpars,windlevel(locs),ka_par_s(locs));

                coeff = (modelpars.maxdensity - snowdensity(locs)) ./ (modelpars.maxdensity - modelpars.density_i);
                kappa_snow = coeff * modelpars.kappa_s_max_density + (1-coeff) * modelpars.kappa_s_density_i;
                g0 = kappa_snow./(state.sdepth(locs)/1000).*(state.T_snowm(locs)+gen.KELVIN-Ts);
                balance1 = Rn_snow0 + g0 + lambdaE0 + H0;

                Ts = Ts0 + 2;
                [~,~,~,~,~,~,~,~,~,H0,lambdaE0,~,~,~,~,Rn_snow0] = get_snow(gen,Ts,T_soil(locs),state.swe(locs),state.liq_swe(locs),state.Q_snow(locs),...
                    state.Q_soil(locs),snowdensity(locs),Rnet_in_g(locs),rain_on_snow_sfc(locs),tsfall(locs),airt(locs),vapp(locs),wind(locs),pres(locs),modelpars,windlevel(locs),ka_par_s(locs));

                coeff = (modelpars.maxdensity - snowdensity(locs)) ./ (modelpars.maxdensity - modelpars.density_i);
                kappa_snow = coeff * modelpars.kappa_s_max_density + (1-coeff) * modelpars.kappa_s_density_i;
                g0 = kappa_snow./(state.sdepth(locs)/1000).*(state.T_snowm(locs)+gen.KELVIN-Ts);
                balance2 = Rn_snow0 + g0 + lambdaE0 + H0;

                slope = (balance2 - balance1)/dT;
                slope(isnan(slope)) = 1;
                balance2(abs(balance1)<thresh) = 0;
                balance1(abs(balance1)>=thresh) = 0;
                balance1 = balance1 + balance2;

                while ~done
                    iter = iter + 1;
                    Ts0 = Ts;

                    dT = -balance1 ./ slope; 
                    dT = max(-2,min(2,dT));
                    Ts = Ts0 + dT;

                    [~,~,~,~,~,~,~,~,~,H0,lambdaE0,~,~,~,~,Rn_snow0] = get_snow(gen,Ts,T_soil(locs),state.swe(locs),state.liq_swe(locs),state.Q_snow(locs),...
                        state.Q_soil(locs),snowdensity(locs),Rnet_in_g(locs),rain_on_snow_sfc(locs),tsfall(locs),airt(locs),vapp(locs),wind(locs),pres(locs),modelpars,windlevel(locs),ka_par_s(locs));

                    coeff = (modelpars.maxdensity - snowdensity(locs)) ./ (modelpars.maxdensity - modelpars.density_i);
                    kappa_snow = coeff * modelpars.kappa_s_max_density + (1-coeff) * modelpars.kappa_s_density_i;
                    g0 = kappa_snow./(state.sdepth(locs)/1000).*(state.T_snowm(locs)+gen.KELVIN-Ts);
                    balance2 = Rn_snow0 + g0 + lambdaE0 + H0;

                    slope = (balance2 - balance1)./dT;
                    slope(isnan(slope)) = 1;
                    balance2(abs(balance1)<thresh) = 0;
                    balance1(abs(balance1)>=thresh) = 0;
                    balance1 = balance1 + balance2;

                    if iter >= niters || max(abs(balance1(:))) < thresh
                        done = 1;
                    end
                end


                coeff = min(1,state.swe(locs) / 50);
                Ts = Ts .* coeff + airt(locs) .* (1-coeff);
                Ts = min(gen.KELVIN,Ts);
                Ts(state.liq_swe(locs) > 0) = gen.KELVIN; 

                [state.swe(locs),state.liq_swe(locs),state.Q_snow(locs),melt(locs),sublimation(locs), ...
                    state.T_snows(locs),state.T_snowm(locs),Q_soil_add(locs),Q_conv(locs),H(locs),lambdaE(locs),ka(locs),...
                    g(locs),Q_m(locs),A_p(locs),Rn_snow(locs)] = get_snow(gen,Ts,T_soil(locs),state.swe(locs),state.liq_swe(locs),...
                    state.Q_snow(locs),state.Q_soil(locs),snowdensity(locs),Rnet_in_g(locs),rain_on_snow_sfc(locs),tsfall(locs),airt(locs),vapp(locs),wind(locs), ...
                    pres(locs),modelpars,windlevel(locs),ka_par_s(locs));

            end

            Temp_g = airt-gen.KELVIN;
            Temp_g_snow = state.T_snowm;
            Temp_g(state.swe>0) = 0;
            Temp_g_snow(state.swe == 0) = 0;
            Temp_g = Temp_g + Temp_g_snow;

            state.depth = ((state.swe+state.liq_swe) ./ snowdensity); state.depth(isnan(state.depth)) = 0;

            %% Energy of Soil Layer

            ice_soil_0 = state.ice_percent_soil .* modelpars.groundlayer_thickness * 1000;

            g_abv = modelpars.kappa_soil ./ (modelpars.groundlayer_thickness / 2) .* (Temp_g - T_soil);
            g_blw = modelpars.kappa_soil ./ (modelpars.dampd-(modelpars.groundlayer_thickness / 2)) .* (modelpars.tempdampd - T_soil);

            g_abv_snow = -(g);
            g_abv(state.swe > 0) = 0;
            g_abv_snow(state.swe <= 0) = 0;
            g_abv = g_abv + g_abv_snow;
            p_dQ = (g_abv + g_blw) * gen.TS;
            gtlocs = state.Q_soil > 0;
            ltlocs = state.Q_soil < 0;
            Q_soil_gtlocs = max(0,state.Q_soil + p_dQ);
            Q_soil_ltlocs = min(0,state.Q_soil + p_dQ);
            Q_soil_gtlocs(gtlocs == 0) = 0; Q_soil_ltlocs(ltlocs == 0) = 0;
            Q_soil = Q_soil_gtlocs + Q_soil_ltlocs;
            residual = state.Q_soil + p_dQ - Q_soil;
            soil_melt = (residual / (modelpars.rhow * modelpars.fusheat)) * (gen.M2MM);

            p_ice_soil = state.smc .* modelpars.groundlayer_thickness * 1000;
            ice_soil_0 = ice_soil_0 - soil_melt;
            ice_soil = max(0,min(p_ice_soil,ice_soil_0));
            residual = ice_soil-ice_soil_0;
            Q_soil = Q_soil + residual * (modelpars.rhow * modelpars.fusheat) / gen.M2MM;
            state.Q_soil = Q_soil;
            state.ice_percent_soil = ice_soil ./ (modelpars.groundlayer_thickness * 1000);

            %% Additional Calculations for output

            Rl_u_g = emiss_g .* modelpars.boltz .* Temp_g.^4;
            Rn_g = Rnet_in_g - Rl_u_g;
            albedo = albedo_g .* (1-vegcover_canopy) + albedo_c .* vegcover_canopy;
            emiss = emiss_g .* (1-vegcover_canopy) + emiss_c .* vegcover_canopy;
            Temp = Temp_g .* (1-vegcover_canopy) + (state.Temp_c - 273.15) .* vegcover_canopy;
            Rn = (Rs_d_direct_top+Rs_d_diffuse) .* (1-albedo) + Rs_scatter + dlwrf - (emiss .* modelpars.boltz .* Temp.^4); 
            airt = airt-gen.KELVIN;
            vapp = vapp/100;
            pres = pres/100;

            % Keep running average of maximum air temperature (for veg. dynamics)
            MaxT = (maxairt-gen.KELVIN);
            state.AvMaxT = (6 * state.AvMaxT + MaxT) / 7;
            T_snows = state.T_snows;
            T_snowm = state.T_snowm;
            T_snows(state.swe == 0) = NaN;
            T_snowm(state.swe == 0) = NaN;
            
            %% Put Data into Output Format

            % Spatial Data (re-expand to matrix form)
            for v = 1:numel(MVars)
                if MVars(v).map
                    var = NaN(sz);
                    eval(['var(~nanlocs) = ' MVars(v).modelvar ';']);
                    eval([MVars(v).modelvar '_expanded = var;']);
                    if strcmp(Info.MapType,'day')
                        eval([MVars(v).Name  '_output = ' MVars(v).Name  '_output + ' MVars(v).modelvar '_expanded * gen.TS/gen.DAY;']);
                    elseif strcmp(Info.MapType,'ts')
                        eval([MVars(v).Name  '_output(:,:,j) = ' MVars(v).modelvar '_expanded;']);
                    end
                end
            end

            % Tabular Data
            if Info.ExportCSV
                for i = 1:numel(SpatialData.DispLocs)
                    if ~isempty(SpatialData.DispLocs(i).pos)
                        for v = 1:numel(MVars)
                            eval(['TabularData(i).' MVars(v).Name '(f,:) = mean(' MVars(v).modelvar '(SpatialData.DispLocs(i).pos(:)));']);
                        end
                    end
                end
            end
        end

        % Only save spatial data of MVars if specified (all MVars will appear
        % as tabular outputs)
        varlist = [];
        for v = 1:numel(MVars)
            if MVars(v).map
                varlist = [varlist ',''' MVars(v).Name  '_output'''];
            end
        end
        
        % Save Restart Files
        if ~exist(restartdirname,'file')
            mkdir(restartdirname)
        end

        % At the end of each day, save model restart file and tabular outputs
        TabularData_disp = TabularData;
        save(restartfilename_next,'state','TabularData_disp');

        % Save the daily output or TS output file
        if ~isempty(varlist)
             eval(['save(modeloutputfile' varlist ');']);
        end

        for v = 1:numel(MVars)
            if MVars(v).map
                if strcmp(Info.MapType,'day')
                    eval([MVars(v).Name  '_output = single(zeros(sz));']);
                elseif strcmp(Info.MapType,'ts')
                    eval([MVars(v).Name  '_output = single(zeros([sz gen.DAY/gen.TS]));']);
                end
            end
        end     
    end
    if Info.OnlyKeepLatestRestart
        dvec = datevec(TS);
        if exist(restartfilename,'file') && ~(dvec(3) == 1)
            delete(restartfilename)
        end
    end
end

if Info.ExportCSV
    if ~exist([Info.DisplayDir filesep 'Tabular'],'file')
        mkdir([Info.DisplayDir filesep 'Tabular']);
    end
    save([Info.DisplayDir filesep 'Tabular' filesep Info.NameIdentifier],'TabularData');
end


function [tsfall_canopy,melt_drip,snow_unload,cansnowstor,acsub] = getcsnowstor(gen,modelpars,Temp_c,Rs_snow,LAI_c,snow,airt,vapp,cansnowstor,wind,rh)
%% Canopy snow storage
% Code to keep track of canopy snow storage
% See Appendix C of Reinhart et al paper
%
% INPUTS
%  state - this is a structure that contains state variables (that is 
%    snowpack water content (swe) and snowpack thermal content (Q)
%  rnet_snow - net radiation at the ground for the current timestep (W/m2)
%  gen.TS - duration of current timestep in seconds
%  metvars - this is a structure containing meteorological variables (wind,
%    airt, snow, rain)
%  hill - this is a structure containing physical terrain characteristics
%  rad - this is a structure containing parameter values
%  veg - this is a structure containing vegetation characteristics
%  modelpars - this is a structure containing parameter values
%  modelpars - this is a structure containing parameter values
%
% OUTPUTS
%  output - this is a variable that can be used for checking
%  state - this module updates the model state and passes it back out
%
% Written by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated 5/2010

cansnowstorcap = 4.4*LAI_c;
L = 0.7 * (cansnowstorcap - cansnowstor) .* (1 - exp(-snow./cansnowstorcap));

tsfall_canopy = snow - L;
cansnowstor = cansnowstor + L;

melt_drip = max(5.8E-5 * (Temp_c - gen.KELVIN) * gen.TS/gen.DAY,0)*1000; % Snow Drip Rate
snow_unload = max(modelpars.snow_unload_par * cansnowstor * gen.TS/gen.DAY,0); % Snow Drip Rate

r = 5E-4;
a = 0.9;
C_e = 0.01 * (cansnowstor./cansnowstorcap) .^ -0.4;               % Calculate Canopy Sublimation
C_e(isnan(C_e)) = 0;
m = modelpars.rhoi * 4/3 * pi * r^3;                        % Mass of Ice Sphere (kg)
rho_v = 0.622*vapp ./ (287 .* Temp_c);    % Water Vapor Density (kg/m3)
S_p = pi * r^2 * (1-a) * Rs_snow;                       % Radiation absorbed by partical (W/m2)
D = 2.06E-5 * (airt/gen.KELVIN) .^ 1.75;                % Diffusivity of air (m2/s)
nu = 1.3E-5;                                            % Kinematic viscosity of air (m2/s)
a_flow = 0.9 * LAI_c;                                   % Canopy flow index
u_c = wind .* exp(-a_flow .* (1-0.6));                  % ventilation velocity
Re = 2 * r * u_c/nu;                                    % Reynolds number
Sh = 1.79 + 0.606 * Re .^ 0.5;                          % Sherwood number
Nu = Sh;                                                % Nusset number
M = 18.01E-3;                                           % Molecular weight of water(kg/mol)
k_t = 0.024;                                             % Thermal conductivity of air (W/m2-K)
R_const = 8314;                                        % Universal gas constant (J/mol-K)

omega = (1 ./ (k_t .* airt .* Nu)) .* ((modelpars.subheat .* M) ./ (R_const .* airt) -1);        % Ventilation factor
dmdt = (2*pi*r*(rh-1) - S_p .* omega) ./ (modelpars.subheat .* omega + 1 ./ (D .* rho_v .* Sh));
psi_s = dmdt/m;   
% Sublimation rate loss coefficient
acsub = C_e .* cansnowstor .* psi_s * gen.TS;
acsub(isnan(acsub)) = 0;
acsub = -acsub;
acsub(acsub>0) = acsub(acsub>0) * modelpars.acsub_mult;

p_canopy_abl = acsub + melt_drip + snow_unload;
deficit = max(1,p_canopy_abl./cansnowstor);
melt_drip = melt_drip./deficit;
snow_unload = snow_unload./deficit;
acsub = acsub./deficit;

cansnowstor = max(0,cansnowstor - acsub - melt_drip - snow_unload); 	% Subtract evaporated snow from the canopy


function [swe,liq_swe,Q,melt,sublimation,Ts,Tm,Q_soil_add,Q_conv,H,lambdaE,ka,g,Q_m,A_p,rnet_snow] = get_snow(gen,Temp_g,T_soil,swe_0,liq_swe_0,Q_0,Q_soil,density,rnet_in_snow,rain,snow,airt,vapp,wind,pres,modelpars,windlevel,ka_par)
%% Run snowmelt program
%
% This is code for a simple snowpack energy balance that allows snow to 
% accumulate (when it is cold) on the land surface until it is warm enough 
% for snowmelt to occur.  It is fully distributed, so that it can compute
% the state of the snowpack across an entire landscape at once.  It is
% designed to work within a fully distributed model with a time loop. In
% addition, IT WAS DEVELOPED TO WORK WITHIN A MODEL RUN AT HOURLY
% TIMESTEPS AND IS UNTESTED AT OTHER TEMPORAL RESOLUTIONS.  Additional
% information can be found in the attached documentation.
%
% Written by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated 5/2010
            
airT_c = airt-gen.KELVIN;
Temp_g = Temp_g-gen.KELVIN;
Ts = Temp_g;

% compute the architectural resistence (bulk transfer coefficient)
ka = ka_par .* wind; 
Tm = (Q_0) ./ (swe_0/1000 * modelpars.specheat_i * modelpars.rhoi);   % Prevent singularities from occurring for extremely shallow snowpacks
coeff = max(0,min(1,(swe_0-10) / 10));
Tm = Tm .* coeff + Ts .* (1-coeff);
Tm(liq_swe_0 > 0) = 0;

% Estimate the richardson number
Ri = (modelpars.g .* (airT_c - Temp_g) .* windlevel) ./ (wind .^2 .* (airT_c + gen.KELVIN));
Ri(Ri>100) = 100;
Ri(Ri<-100) = -100;

% Adjust the bulk transfer coefficient
tmp = ka ./ (1 + 10 * Ri);
tmp2 = ka .* (1 - 10 * Ri);
tmp(Ri < 0) = 0;
tmp2(Ri >= 0) = 0;
ka = tmp + tmp2;
ka(isnan(ka)) = 0;

% Compute the snow surface vapor pressure (assume to be saturated at snow temperature)
svapp = 0.6108*exp(17.27*Ts./(237.3+Ts))*1000; % in Pa

% Compute the latent heat flux
lambdaE = modelpars.rhoa * modelpars.subheat .* ka .* (0.622*(vapp-svapp)./pres) * modelpars.EMult;
lambdaE(isnan(lambdaE)) = 0;

sublimation = -(lambdaE / modelpars.subheat) * (gen.M2MM/modelpars.rhoi) * gen.TS;
evaporation = -(lambdaE / modelpars.evapheat) * (gen.M2MM/modelpars.rhoi) * gen.TS;
sublimation(liq_swe_0 > 0) = 0;
evaporation(liq_swe_0 == 0) = 0;
sublimation = sublimation + evaporation;
swe_1 = max(0,swe_0 - sublimation);
sublimation = swe_0-swe_1;
lambdaE = -sublimation * modelpars.subheat * modelpars.rhoi / (gen.M2MM * gen.TS);

% Compute the sensible heat flux
Cx = 2;     % J-kg/m2-s
H = (modelpars.rhoa * modelpars.specheat_a .* ka + Cx) .* (airT_c - Ts) * modelpars.HMult;
H(isnan(H)) = 0;

% Compute snow surface net radiation

rnet_snow = rnet_in_snow - modelpars.emiss_snow .* modelpars.boltz .* (Ts+gen.KELVIN).^4;

% Compute heat from precip

A_p = modelpars.specheat_w * modelpars.rhow * ((rain+snow) / (gen.M2MM * gen.TS) ) .* (airT_c - Ts);

% Compute the ground heat flux

coeff = (modelpars.maxdensity - density) ./ (modelpars.maxdensity - modelpars.density_i);
kappa_snow = coeff * modelpars.kappa_s_max_density + (1-coeff) * modelpars.kappa_s_density_i;
g = (kappa_snow ./ (((swe_0/1000)./density)/2 + modelpars.groundlayer_thickness/2)) .* (T_soil-Tm);          % McCumber and Pielke, 1981

dQ = (rnet_snow + H + lambdaE + A_p + g) * gen.TS;


% sublimation = sublimation./deficit;
swe_0 = max(0,swe_0 - sublimation);

% Update the snowpack internal energy (use logical indexing for distributed model)
Q_m = Q_0 + dQ;
Q_m(Q_m<0) = 0;
Q = Q_0 + dQ - Q_m;

% Compute snowmelt and update snow water equivalent (use logical indexing for distributed model)
melt = (Q_m / (modelpars.rhow * modelpars.fusheat)) * (gen.M2MM);   % in mm/timestep
swe = max(0, swe_0 - melt);
melt = swe_0 - swe; melt(melt < 0) = 0;

liq_swe_cap = swe .* modelpars.liq_swe_cap_mult;
liq_swe = liq_swe_0 + melt + rain;
melt = max(0,liq_swe - liq_swe_cap);
liq_swe = liq_swe - melt;
% Recalculate melt heat so it represents melt coming out of pack
Q_m_new = melt .* (modelpars.rhow * modelpars.fusheat) / (gen.M2MM * gen.TS);
Q_conv_m = (Q_m - Q_m_new) / gen.TS;
Q_m = Q_m_new;

% A_pw_pack = modelpars.specheat_w * modelpars.rhow * (liq_swe / (gen.M2MM * gen.TS) ) .* (0 - Tm);
p_ref_pack = -Q ./ (modelpars.rhow * modelpars.fusheat) * (gen.M2MM);
a_ref_pack = min(p_ref_pack,liq_swe);
liq_swe = liq_swe - a_ref_pack;
swe = max(0,swe + a_ref_pack);
Q_new = Q + (modelpars.rhow * modelpars.fusheat * a_ref_pack / (gen.M2MM));
Q_conv = (Q - Q_new) / gen.TS;
Q = Q_new;

% Find temperature of the snowpack (based on new cold content)
Q(swe == 0) = 0;

Tm = (Q) ./ (swe/1000 * modelpars.specheat_i * modelpars.rhoi); % Prevent singularities from occurring for extremely shallow snowpacks
coeff = max(0,min(1,(swe-10) / 10));
Tm = Tm .* coeff + Ts .* (1-coeff);

Tm(liq_swe > 0) = 0;
Tm(isnan(Tm)) = 0;
% Recalculate cold content and put remainder in the ground
Q_new = Tm .* (swe/1000 * modelpars.specheat_i * modelpars.rhoi);

Q_soil_add = (Q_new - (Q_0 + (rnet_snow + H + lambdaE + A_p + g - Q_conv_m - Q_conv) * gen.TS)) / gen.TS;
Q_conv = Q_conv_m + Q_conv;
Q = max(Q_new,-modelpars.specheat_i*modelpars.rhow*swe./gen.M2MM*20);
