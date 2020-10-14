function [prec,airt,pres,vapp,ugrd,vgrd,dswrf,dlwrf] = get_forcing_data(Info,SpatialData,TS,VarInclude)
% Function to load the forcing data and interpolate to the model grid for a
% single timestep
%
% Inputs: Info - Structure containing information about the model run
%         SpatialData - Structure containing spatial data used in the model
%         TS - Timestep for the day to retrieve forcing data for
%         VarInclude - Variables to retrieve (right now, precip, airtemp,
%           pressure, vapor pressure, u wind, v wind, downward shortwave and
%           longwave radiation
% Outputs: Maps of the forcing data, interpolated to the model grid
%
% Timestep information
dvec = datevec(TS);
yr = num2str(dvec(1));
mo = num2str(dvec(2)); if numel(mo) < 2, mo = ['0' mo]; end
dy = num2str(dvec(3)); if numel(dy) < 2, dy = ['0' dy]; end
hr = num2str(dvec(4)); if numel(hr) < 2, hr = ['0' hr]; end
mn = num2str(dvec(5)); if numel(mn) < 2, mn = ['0' mn]; end

% Unpack spatial data
elev = SpatialData.z;
X = SpatialData.LonMap;
Y = SpatialData.LatMap;

vapp = NaN(size(elev));

% Loop through all of the variables
for i = 1:numel(Info.Forcing.var)
    % Different variables can come from different sources (see ProgramPars)
    ForcingDir = Info.ForcingDirs{Info.Forcing.source{i}};
    % They can have different grids
    lon2d = Info.Forcing.lon2d{Info.Forcing.source{i}};
    lat2d = Info.Forcing.lat2d{Info.Forcing.source{i}};
    Z = Info.Forcing.Z{Info.Forcing.source{i}};
    sz = Info.Forcing.sz{Info.Forcing.source{i}};
    % And be in different places in the NetCDF files
    start = Info.Forcing.start{Info.Forcing.source{i}};
    count = Info.Forcing.count{Info.Forcing.source{i}};
    stride = Info.Forcing.stride{Info.Forcing.source{i}};
    % Interpolation methods can be different
    useLapseRate = Info.Forcing.useLapseRate{i};
    interpMethod = Info.Forcing.interpMethod{i};
    % And they can come from different variables in the files
    var = Info.Forcing.var{i};
    fprefix = Info.Forcing.fprefix{i};
   
    if VarInclude(i) == 1
        % Read the forcing information
        data = ncread([ForcingDir filesep yr filesep mo filesep dy filesep fprefix '.' yr mo dy hr mn '.nc'],var,start,count,stride);

        % If specified, interpolate based on local lapse rates instead of
        % doing basic interpolation
        if useLapseRate
            warning off all
            M= NaN(size(data));
            B = NaN(size(data));
            for r = 2:sz(1)-1
                for c = 2:sz(2)-1
                    y = single([data(r-1,c+1); data(r-1,c); data(r-1,c-1); data(r,c+1); data(r,c); data(r,c-1); data(r+1,c+1); data(r+1,c); data(r+1,c-1)]);
                    x = single([Z(r-1,c+1); Z(r-1,c); Z(r-1,c-1); Z(r,c+1); Z(r,c); Z(r,c-1); Z(r+1,c+1); Z(r+1,c); Z(r+1,c-1)]);
                    if sum(isnan(y)) < 5
                        nanlocs = isnan(x) | isnan(y);
                        MB = [ones(length(x(~nanlocs)),1) x(~nanlocs)] \ y(~nanlocs);
                        M(r,c) = MB(2);
                        B(r,c) = MB(1);
                    end
                end
            end
            M(isnan(data)) = NaN;
            B(isnan(data)) = NaN;
            M(data == 0 | isnan(M)) = 0;
            B(data == 0 | isnan(B)) = 0;
            warning on all
            % Interpolate both m and b
            M_grid = griddedInterpolant(single(lon2d),single(lat2d),single(M),interpMethod);
            B_grid = griddedInterpolant(single(lon2d),single(lat2d),single(B),interpMethod);
            newdata = M_grid(X,Y).*single(elev) + B_grid(X,Y);
        else
            % Else, just perform basic interpolation on the data
            data_grid = griddedInterpolant(single(lon2d),single(lat2d),single(data),interpMethod);
            newdata = data_grid(X,Y);
        end

        % Fill the output variables with the forcing data
        switch var
            case 'PRES'
                pres = newdata;
            case 'TMP'
                airt = newdata;
            case 'UGRD'
                ugrd = newdata;
            case 'VGRD'
                vgrd = newdata;
            case 'SPFH'
                spfh = newdata;
            case 'APCP'
                prec = newdata;
            case 'DLWRF'
                dlwrf = newdata;
            case 'DSWRF'
                dswrf = newdata;
        end
       
        if strcmp(var,'SPFH')
            vapp = pres.*spfh ./ (0.622+spfh);
        end
    else   
        switch var
            case 'PRES'
                pres = NaN(size(elev));
            case 'TMP'
                airt = NaN(size(elev));
            case 'UGRD'
                ugrd = NaN(size(elev));
            case 'VGRD'
                vgrd = NaN(size(elev));
            case 'SPFH'
                spfh = NaN(size(elev));
            case 'APCP'
                prec = NaN(size(elev));
            case 'DLWRF'
                dlwrf = NaN(size(elev));
            case 'DSWRF'
                dswrf = NaN(size(elev));
        end
    end
end

