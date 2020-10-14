function [R_tot,Id_tot,cos_i] = solarradiation_direct(slop,asp,lat,Doy,tstart,tend,n,Alphas,AlphaMat)
% PUPROSE: Calculate incoming solar radiation for SnowPALM
% -------------------------------------------------------------------
% Note: This function is based on the file exchange submission 'solarrradiation' by Felix Hebeler. 
% It has been modified for SnowPALM

% USAGE: [R_tot,Id_tot,cos_i] = solarradiation_direct(slop,asp,lat,Doy,tstart,tend,n,Alphas,AlphaMat)
% where: slop is the slope map
%        asp is the aspect map
%        lat is the map of latitudes
%        Doy is the day of year (scaler)
%        tstart is the starting time for the calculation
%        tend is the ending time for the calculation
%        n is the timestep of calculation: 1=hourly, 0.5=30min, 2=2hours etc
%        Alphas is a vector giving the directions that are in AlphaMat
%        AlphaMat: 3-D matrix of horizon angles for each direction 
%        cs is the cellsize in meters
%        r is the ground reflectance (global value or map, default is 0.2)
%        R_tot is a map of direct solar radiation (W/m2)
%        Id_tot is a map of diffuse solar radiation (W/m2)
%        cos_i is the solar incedence angle
%
% See also: GRADIENT, CART2POL
%
% Note: Follows the approach of Kumar et al 1997. Calculates clear sky
%       radiation corrected for the incident angle (selfshading) plus
%       diffuse and reflected radiation. Insolation is depending on time of year (and day), 
%       latitude, elevation, slope and aspect. 
%       Relief shading is not considered.
%       Script uses simple unweighed gradient of 4 nearest neighbours for slope
%       calculation.
%
% Reference: Kumar, L, Skidmore AK and Knowles E 1997: Modelling topographic variation in solar radiation in 
%            a GIS environment. Int.J.Geogr.Info.Sys. 11(5), 475-497
%
% Felix Hebeler, Dept. of Geography, University Zurich, May 2008.
% Modified by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated July 2017


%% parameters
%It ;               % total hours of daily sunshine (calculated inline)
%M ;                % air mass ratio parameter (calculated inline)
%r = 0.20;          % ground reflectance coefficient (more sensible to give as input)
%L=lat;             %latitude
% % PDB: Modified to accept n as an input instead of being defined here
% n = 1;              % timestep of calculation over sunshine hours: 1=hourly, 0.5=30min, 2=2hours etc
tau_a    = 365;     %length of the year in days
S0 = 1367;          % solar constant W m^-2   default 1367

dr= 0.0174532925;   % degree to radians conversion factor

%%  convert factors
% % PDB: slop and asp are given as inputs instead of calculated here
% [slop,asp]=get_ders(dem,cs);   % calculate slope and aspect in radians using given cellsize cs
% % PDB: Modified to accept lat as an input
% [dummy,L]=meshgrid(1:size(dem,2),lat);   % grid latitude
% clear dummy;
L=lat;
L=L*dr;                     % convert to radians
fcirc = 360*dr; % 360 degrees in radians

%% some setup calculations
% % PDB: No longer outputting srad, but instead outputting its components
% srad=0;
sinL=sin(L);
cosL=cos(L);
tanL=tan(L);
sinSlop=sin(slop);
cosSlop=cos(slop);
cosSlop2=cosSlop.*cosSlop;
sinAsp=sin(asp);
cosAsp=cos(asp);
term1 = ( sinL.*cosSlop - cosL.*sinSlop.*cosAsp);
term2 = ( cosL.*cosSlop + sinL.*sinSlop.*cosAsp);
term3 = sinSlop.*sinAsp;

% % PDB: Do not loop through the year, accept a day of year instead
% for d = 1:tau_a; 
    d = Doy;
    %display(['Calculating melt for day ',num2str(d)])  
    % clear sky solar radiation
    I0 = S0 * (1 + 0.0344*cos(fcirc*d/tau_a)); % extraterr rad per day     
    % sun declination dS
    dS = 23.45 * dr* sin(fcirc * ( (284+d)/tau_a ) ); %in radians, correct/verified
    % % PDB: Revised method to get hour angle
    % % angle at sunrise/sunset
    % % t = 1:It; % sun hour    
    % hsr = real(acos(-tanL*tan(dS)));  % angle at sunrise
    % % this only works for latitudes up to 66.5 deg N! Workaround:
    % % hsr(hsr<-1)=acos(-1);
    % % hsr(hsr>1)=acos(1);
    % It=round(12*(1+mean(hsr(:))/pi)-12*(1-mean(hsr(:))/pi)); % calc daylength
    %%  daily loop
    % % PDB: Instead of looping over sunshine hours, we go from t_start to
    % t_end
%     I=0;
%     for t=1:n:It % loop over sunshine hours
    R_tot = 0;
    Rdif_tot = 0;
    c = 0;
    for t=tstart:n:tend
        c = c+1;
        % if accounting for shading should be included, calc hillshade here
        % hourangle of sun hs 
        % % PDB: Revised method to get hour angle
        %  hs=hsr-(pi*t/It);               % hs(t)
        hs = (pi/12)*((t-0.5)-12);                      	% Hour Angle (PDB)
        %solar angle and azimuth
        %alpha = asin(sinL*sin(dS)+cosL*cos(dS)*cos(hs));% solar altitude angle
        sinAlpha = sinL.*sin(dS)+cosL.*cos(dS).*cos(hs);
        %alpha_s = asin(cos(dS)*sin(hs)/cos(alpha)); % solar azimuth angle
        % correction  using atmospheric transmissivity taub_b
        M=sqrt(1229+((614.*sinAlpha)).^2)-614.*sinAlpha; % Air mass ratio
        tau_b = 0.56 * (exp(-0.65*M) + exp(-0.095*M));
        tau_d = 0.271-0.294*tau_b; % radiation diffusion coefficient for diffuse insolation
        tau_r = 0.271+0.706*tau_b; % reflectance transmitivity
        % correct for local incident angle
        cos_i = (sin(dS).*term1) + (cos(dS).*cos(hs).*term2) + (cos(dS).*term3.*sin(hs));
        Is = I0 * tau_b; % potential incoming shortwave radiation at surface normal (equator)
        % R = potential clear sky solar radiation W m2
        R = Is .* cos_i;
        R(R<0)=0;  % kick out negative values
        Id = I0 .* tau_d .* cosSlop2./ 2.*sinAlpha; %diffuse radiation;
        % % PDB: Do not compute reflectance here, output components separately
%         Ir = I0 .* r .* tau_r .* sinSlop2./ 2.* sinAlpha; % reflectance
%         R= R + Id + Ir;
%         R(R<0)=0; 
%         I=I+R;% solar radiation per day (sunshine hours)  
%      end % end of sun hours in day loop
% %%  add up radiation part melt for every day
%     srad = srad + I;
% end   % end of days in year loop
    
    

    %% Modified by Patrick Broxton (10/2010)
    
    % Get rid of bad values of Id
    Id(Id < 0) = 0;     % kick out negative values
    Id = Id ./ (R+Id);
    Id(Id <= 0) = 1;
    Id(isnan(Id)) = 1;  % Clear sky diffuse fraction
    Rdif = R .* Id; % Clear sky diffuse fraction
    % Calculate shadowing effects (sun below horizon)
    alpha = asin(sinAlpha);
    if ~isempty(Alphas)
        alpha_s = -asin(cos(dS).*sin(hs)./cos(alpha));            % Solar Azimuth Angle
        alpha_sx = -(cos(hs).*cos(dS).*sinL)+(cosL.*sin(dS));      % Solar Declination
        alpha_s2 = pi-alpha_s;
        alpha_s(alpha_sx>0) = 0;
        alpha_s2(alpha_sx<0) = 0;
        alpha_s = alpha_s + alpha_s2;
        alpha_s2 = alpha_s + 2*pi;
        alpha_s(alpha_s<0) = 0;
        alpha_s2(alpha_s>0) = 0;
        alpha_s = alpha_s + alpha_s2; 
        dirs = find(alpha_s(1) > Alphas & alpha_s(1) <= Alphas+1,1,'last');
        min_alphas = zeros(size(alpha_s));
        if ndims(AlphaMat) == 3
            min_alpha = double(AlphaMat(:,:,dirs))/(255/(pi/2));
        else
            min_alpha = double(AlphaMat(:,dirs))/(255/(pi/2));
        end
        min_alpha(dirs ~= dirs) = 0;
        min_alphas = min_alphas + min_alpha;
        R(alpha<min_alphas) = 0;
    end
    R_tot = R_tot + R;
    Rdif_tot = Rdif_tot + Rdif;
end
R_tot = R_tot/c;
Rdif_tot = Rdif_tot/c;
Id_tot = Rdif_tot ./ R_tot;
%% End Modification

% %%
% function [grad,asp] = get_ders(dem,cs)
% % calculate slope and aspect (deg) using GRADIENT function
% [fx,fy] = gradient(dem,cs,cs); % uses simple, unweighted gradient of immediate neighbours
% [asp,grad]=cart2pol(fy,fx); % convert to carthesian coordinates
% grad=atan(grad); %steepest slope
% asp=asp.*-1+pi; % convert asp 0 facing south
