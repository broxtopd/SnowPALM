function [ForcingVals] = get_daily_nldas_forcing(Info,SpatialData,TS,VarInclude)
% Reads gridded forcing data, and interpolates it to the model grid, and
% returns a structure containing daily forcing data
%
% Inputs: Info - Structure containing information about the model run
%         SpatialData - Structure containing spatial data used in the model
%         TS - Timestep for the day to retrieve forcing data for
%         VarInclude - Variables to retrieve (right now, precip, airtemp,
%           pressure, vapor pressure, u wind, v wind, downward shortwave and
%           longwave radiation
% Outputs: ForcingVals - Structure containing forcing data for a given day
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated June 2016

% Timestep information
dvec = datevec(TS);
yr = num2str(dvec(1));
mo = num2str(dvec(2)); if numel(mo) < 2, mo = ['0' mo]; end
dy = num2str(dvec(3)); if numel(dy) < 2, dy = ['0' dy]; end
% Forcing Directory
fdir = [Info.ModelStatesDir filesep 'Forcing Files' filesep Info.NameIdentifier filesep yr filesep mo];
fname = [dy '.mat'];

% If the forcing intepolation file is not saved for a particular day ...
if ~exist([fdir filesep fname],'file')
    c = 0;
    for hour = 0:23
        c = c+1;
        ts = TS + hour/24;
        % Function to load the forcing data and interpolate to the model grid
        [prec,airt,pres,vapp,ugrd,vgrd,dswrf,dlwrf] = get_forcing_data(Info,SpatialData,ts,VarInclude);
        PREC(:,c) = prec;
        AIRT(:,c) = airt;
        PRES(:,c) = pres;
        VAPP(:,c) = vapp;
        UGRD(:,c) = ugrd;
        VGRD(:,c) = vgrd;
        DSWRF(:,c) = dswrf;
        DLWRF(:,c) = dlwrf;
    end

    % If specified, save the forcing interpolation file (so we don't need
    % to interpolate again
    if Info.SaveForcingInterpFiles
        if ~exist(fdir,'file')
            mkdir(fdir);
        end
        save([fdir filesep fname],'PREC','AIRT','PRES','VAPP','UGRD','VGRD','DSWRF','DLWRF');
    end 
% else, load the saved file
else
	load([fdir filesep fname]); 
end
            
% Put the forcing data into a structure that is used by SnowPALM
c = 0;
for hour = 0:23
    c = c+1;
    ForcingVals(c).prec = PREC(:,c);
    ForcingVals(c).airt = AIRT(:,c);
    ForcingVals(c).pres = PRES(:,c);
    ForcingVals(c).vapp = VAPP(:,c);
    ForcingVals(c).ugrd = UGRD(:,c);
    ForcingVals(c).vgrd = VGRD(:,c);
    ForcingVals(c).dswrf = DSWRF(:,c);
    ForcingVals(c).dlwrf = DLWRF(:,c);
end


