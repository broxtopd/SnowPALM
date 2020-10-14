function get_prism_corrected_nldas(Domain,StartMonth,StartYear,EndMonth,EndYear)
% Download and process PRISM-corrected NLDAS forcing data.  Output data
% will
% Usage: get_prism_corrected_nldas(Domain,StartMonth,StartYear,EndMonth,EndYear)
% Inputs: - Domain (a region defined in domains.m for which to download
%               forcing data for)
%           StartMonth: month of starting date
%           StartYear: year of starting date
%           EndMonth: month of ending date
%           EndYear: year of ending date
% Example
% get_prism_corrected_nldas('CO_NM',1,2013,12,2013)

NLDASDataLoc = 'https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_FORA0125_H.002';
NLDASUsername = 'YOURUSERNAME';
NLDASPassword = 'YOURPASSWORD';
PRISMDataLoc = 'ftp://prism.nacse.org/monthly';

prism_ppt_version = '3';            % Version of PRISM Precipitation data                   
prism_tmean_version = '2';          % Version of PRISM Temperature data

name = ['prism_corrected_nldas_' Domain];               % Name of ouput forcing directory
[ulx,uly,lrx,lry,outdir,demfile] = domains(Domain);     % Get information about the domains

% Timestamps of start and end dates
STStamp = datenum(StartYear,StartMonth,1);
ETStamp = datenum(EndYear,EndMonth,1);
% Forcing bounds
n_ulx = str2double(ulx);
n_uly = str2double(uly);
n_lrx = str2double(lrx);
n_lry = str2double(lry);
% Days in month
DIM_i = [31 28 31 30 31 30 31 31 30 31 30 31];

% Download data to a temporary directory
tmp_dir = tempname;
mkdir(tmp_dir);
if ~exist([outdir filesep name],'file')
    mkdir([outdir filesep name]);
end

% Loop through desired dates
d = 0;
for TS = STStamp:ETStamp
    dvec = datevec(TS);
    d = d+1;
    
    % On the first day of the month, download a bunch of NLDAS data and
    % some PRISM data
    if dvec(3) == 1
        % Get strings for year, month, day, and day of year
        doy = TS - datenum(dvec(1),1,1) + 1;
        yr = num2str(dvec(1));
        mo = num2str(dvec(2)); if numel(mo) < 2, mo = ['0' mo]; end
        dy = num2str(dvec(3)); if numel(dy) < 2, dy = ['0' dy]; end
        doy_str = num2str(doy);
        if numel(doy_str) < 3, doy_str = ['0' doy_str]; end
        if numel(doy_str) < 3, doy_str = ['0' doy_str]; end
        
        if d == 1
            % On the first iteration, download some dummy data, get latitude,
            %  longitude information, also output resized DEM data
            hr = '00';
            url = [NLDASDataLoc '/' yr '/' doy_str '/NLDAS_FORA0125_H.A' yr mo dy '.' hr '00.002.grb'];
            of = [tmp_dir '/' hr '.grb'];
            % Download NLDAS data
            eval(['!wget --user ' NLDASUsername ' --password ' NLDASPassword ' -O "' of '" "' url '"']);
            % Convert to NetCDF so we can easily read the file
            eval(['!gdal_translate -of NetCDF "' of '" "' tmp_dir '/' hr '.nc"']);
            lat = ncread([tmp_dir '/' hr '.nc'],'lat');
            lon = ncread([tmp_dir '/' hr '.nc'],'lon');
            locx = lon > n_ulx & lon < n_lrx;
            locy = lat > n_lry & lat < n_uly; 
            sz = [numel(lon(locx)) numel(lat(locy))];

            % Cut out and resize the DEM to the correct dimensions
            eval(['!gdal_translate -outsize ' num2str(sz(1)) ' '  num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' demfile '" "' tmp_dir '/dem_sub.tif"']);
            dem_small = imread([tmp_dir '/dem_sub.tif']);
            Z = fliplr(dem_small');
            
            % Write DEM, Latitude, Longitude data to NetCDF file
            ncfile = [outdir filesep name filesep 'dem.nc'];
            if exist(ncfile,'file')
                delete(ncfile);
            end
                
            nccreate(ncfile,'Elev','Dimensions', {'lon', sz(1),'lat', sz(2)},'DataType','single');
            ncwriteatt(ncfile,'/Elev','units','meters');
            ncwrite(ncfile,'Elev',Z);
                
            nccreate(ncfile,'lat','Dimensions', {'lat', sz(2)},'DataType','single');
            ncwriteatt(ncfile,'/lat','units','degrees_north');
            ncwrite(ncfile,'lat',lat(locy));

            nccreate(ncfile,'lon','Dimensions', {'lon', sz(1)},'DataType','single');
            ncwriteatt(ncfile,'/lon','units','degrees_east');
            ncwrite(ncfile,'lon',lon(locx));
        end

        disp(['Getting PRISM data for ' datestr(TS,'mmm-yyyy')]);
     
        % Download the PRISM data
        try
            % First try downloading the stable data (available further 
                % back in time than ~6 months in the past)
            of = [tmp_dir '/PRISM_ppt_stable_4kmM' prism_ppt_version '_' yr mo '_bil.zip'];
            url = [PRISMDataLoc '/ppt/' yr '/PRISM_ppt_stable_4kmM' prism_ppt_version '_' yr mo '_bil.zip'];
            eval(['!wget -O "' of '" "' url '"']);
            unzip(of,tmp_dir);
            eval(['!gdal_translate -outsize ' num2str(sz(1)) ' ' num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' tmp_dir '/PRISM_ppt_stable_4kmM' prism_ppt_version '_' yr mo '_bil.bil" "' tmp_dir '/PRISM_ppt_stable_4kmM' prism_ppt_version '_' yr mo '_bil.tif"']);
            ppt = imread([tmp_dir '/PRISM_ppt_stable_4kmM' prism_ppt_version '_' yr mo '_bil.tif']);
            ppt(ppt == -9999) = NaN;
            ppt = fliplr(ppt');


            of = [tmp_dir '/PRISM_tmean_stable_4kmM' prism_tmean_version '_' yr mo '_bil.zip'];
            url = [PRISMDataLoc '/tmean/' yr '/PRISM_tmean_stable_4kmM' prism_tmean_version '_' yr mo '_bil.zip'];
            eval(['!wget -O "' of '" "' url '"']);
            unzip(of,tmp_dir);
            unzip([tmp_dir '/PRISM_tmean_stable_4kmM' prism_tmean_version '_' yr mo '_bil.zip'],tmp_dir);
            eval(['!gdal_translate -outsize ' num2str(sz(1)) ' ' num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' tmp_dir '/PRISM_tmean_stable_4kmM' prism_tmean_version '_' yr mo '_bil.bil" "' tmp_dir '/PRISM_tmean_stable_4kmM' prism_tmean_version '_' yr mo '_bil.tif"']);
            tmean = imread([tmp_dir '/PRISM_tmean_stable_4kmM' prism_tmean_version '_' yr mo '_bil.tif']);
            tmean(tmean == -9999) = NaN;
            tmean = fliplr(tmean');
        catch
            try 
                % Then try downloading the provisional data (available further 
                % back in time than ~1 months in the past)
                of = [tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.zip'];
                url = [PRISMDataLoc '/ppt/' yr '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.zip'];
                eval(['!wget -O "' of '" "' url '"']);
                unzip(of,tmp_dir);
                eval(['!gdal_translate -outsize ' num2str(sz(1)) ' ' num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.bil" "' tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.tif"']);
                ppt = imread([tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.tif']);
                ppt(ppt == -9999) = NaN;
                ppt = fliplr(ppt');


                of = [tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.zip'];
                url = [PRISMDataLoc '/tmean/' yr '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.zip'];
                eval(['!wget -O "' of '" "' url '"']);
                unzip(of,tmp_dir);
                unzip([tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.zip'],tmp_dir);
                eval(['!gdal_translate -outsize ' num2str(sz(1)) ' ' num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.bil" "' tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.tif"']);
                tmean = imread([tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.tif']);
                tmean(tmean == -9999) = NaN;
                tmean = fliplr(tmean');
            catch
                % Then try downloading the early data (more recent than ~1 month in the past)
                of = [tmp_dir '/PRISM_ppt_early_4kmM' prism_ppt_version '_' yr mo '_bil.zip'];
                url = [PRISMDataLoc '/ppt/' yr '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.zip'];
                eval(['!wget -O "' of '" "' url '"']);
                unzip(of,tmp_dir);
                eval(['!gdal_translate -outsize ' num2str(sz(1)) ' ' num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.bil" "' tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.tif"']);
                ppt = imread([tmp_dir '/PRISM_ppt_provisional_4kmM' prism_ppt_version '_' yr mo '_bil.tif']);
                ppt(ppt == -9999) = NaN;
                ppt = fliplr(ppt');


                of = [tmp_dir '/PRISM_tmean_early_4kmM' prism_tmean_version '_' yr mo '_bil.zip'];
                url = [PRISMDataLoc '/tmean/' yr '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.zip'];
                eval(['!wget -O "' of '" "' url '"']);
                unzip(of,tmp_dir);
                unzip([tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.zip'],tmp_dir);
                eval(['!gdal_translate -outsize ' num2str(sz(1)) ' ' num2str(sz(2)) ' -projwin ' ulx ' ' uly ' ' lrx ' ' lry ' "' tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.bil" "' tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.tif"']);
                tmean = imread([tmp_dir '/PRISM_tmean_provisional_4kmM' prism_tmean_version '_' yr mo '_bil.tif']);
                tmean(tmean == -9999) = NaN;
                tmean = fliplr(tmean');
            end
        end

        % Number of days in the month
        DIM = DIM_i(dvec(2));
        if dvec(2) == 2 && (mod(dvec(1),4) == 0 && (mod(dvec(1),100) > 0 || mod(dvec(1),400) == 0))
            DIM = DIM + 1;
        end

        % Download the NLDAS data
        dvec_i = dvec;
        airt = zeros([sz(1) sz(2) 24*DIM]);
        spfh = zeros([sz(1) sz(2) 24*DIM]);
        pres = zeros([sz(1) sz(2) 24*DIM]);
        uwind = zeros([sz(1) sz(2) 24*DIM]);
        vwind = zeros([sz(1) sz(2) 24*DIM]);
        dlwrf = zeros([sz(1) sz(2) 24*DIM]);
        apcp = zeros([sz(1) sz(2) 24*DIM]);
        dswrf = zeros([sz(1) sz(2) 24*DIM]);
        airt_min = zeros([sz(1) sz(2) DIM]);
        airt_max = zeros([sz(1) sz(2) DIM]);
        c = 0; e = 0;
        for dim = 1:DIM
            e = e+1;
            if dim > 1
                dvec(3) = dvec(3) + 1;
                doy = doy + 1;
            end
            dy = num2str(dvec(3)); if numel(dy) < 2, dy = ['0' dy]; end
            doy_str = num2str(doy);
            if numel(doy_str) < 3, doy_str = ['0' doy_str]; end
            if numel(doy_str) < 3, doy_str = ['0' doy_str]; end
                
            disp(['Getting NLDAS data for ' dy '-' datestr(TS,'mmm-yyyy')]);
            for hour = 0:23
                c = c+1;
                hr = num2str(hour);
                if numel(hr) == 1, hr = ['0' hr]; end
                   
                % Try up to 5 times to get the NLDAS data (because it 
                % sometimes fails), read into 3d matrix (all hours in
                % month)
                for try_iter = 1:5
                    try
                        if ~exist([hr '.nc'],'file')
                            
                            url = [NLDASDataLoc '/' yr '/' doy_str '/NLDAS_FORA0125_H.A' yr mo dy '.' hr '00.002.grb'];
                            of = [tmp_dir '/' hr '.grb'];
                            eval(['!wget --user ' NLDASUsername ' --password ' NLDASPassword ' -O "' of '" "' url '"']);
                            eval(['!gdal_translate -of NetCDF "' of '" "' tmp_dir '/' hr '.nc"']);
                        end
                        
                        AIRT = ncread([tmp_dir '/' hr '.nc'],'Band1');       airt(:,:,c) = AIRT(locx,locy);
                        SPFH = ncread([tmp_dir '/' hr '.nc'],'Band2');       spfh(:,:,c) = SPFH(locx,locy);
                        PRES = ncread([tmp_dir '/' hr '.nc'],'Band3');       pres(:,:,c) = PRES(locx,locy);
                        UWIND = ncread([tmp_dir '/' hr '.nc'],'Band4');      uwind(:,:,c) = UWIND(locx,locy);
                        VWIND = ncread([tmp_dir '/' hr '.nc'],'Band5');      vwind(:,:,c) = VWIND(locx,locy);
                        DLWRF = ncread([tmp_dir '/' hr '.nc'],'Band6');      dlwrf(:,:,c) = DLWRF(locx,locy);
                        APCP = ncread([tmp_dir '/' hr '.nc'],'Band10');      apcp(:,:,c) = APCP(locx,locy);
                        DSWRF = ncread([tmp_dir '/' hr '.nc'],'Band11');     dswrf(:,:,c) = DSWRF(locx,locy);
                        % If success, then kick out of loop
                        break
                    catch
                        if try_iter < 5
                            disp('Problem Downloading file, retrying...');
                        else
                            disp('Problem Downloading file, rerun to retry');
                        end
                    end
                end
            end
            airt_min(:,:,e) = min(airt(:,:,c-23:c),[],3);
            airt_max(:,:,e) = max(airt(:,:,c-23:c),[],3);
            delete([tmp_dir '/*']);
        end
        airt(airt == 9999) = NaN;
        spfh(spfh == 9999) = NaN;
        pres(pres == 9999) = NaN;
        uwind(uwind == 9999) = NaN;
        vwind(vwind == 9999) = NaN;
        dlwrf(dlwrf == 9999) = NaN;
        apcp(apcp == 9999) = NaN;
        dswrf(dswrf == 9999) = NaN;
        airt_min(airt_min >= 9999) = NaN;
        airt_max(airt_max >= 9999) = NaN; 
        
        % Process the monthly data, compare the PRISM data to the NLDAS
        % data, apply a bias correction to make monthly NLDAS data the same
        % as PRISM (match cumulative precip, monthly average temp)
        
        ppt_bias = max(0,min(10,sum(apcp,3) ./ ppt));
        ppt_bias(isnan(ppt)) = NaN;
        ppt_bias = repmat(ppt_bias,[1 1 24*DIM]);
        apcp_corr = apcp ./ ppt_bias;
        
        tmean_bias = max(-10,min(10,mean((airt_max+airt_min)/2,3) - tmean));
        tmean_bias(isnan(tmean)) = NaN;
        tmean_bias = repmat(tmean_bias,[1 1 24*DIM]);
        airt_corr = airt + tmean_bias;

        % Write out data for all days in the month
        dvec = dvec_i;
        c = 0;
        for dim = 1:DIM
            if dim > 1
                dvec(3) = dvec(3) + 1;
            end
            dy = num2str(dvec(3)); if numel(dy) < 2, dy = ['0' dy]; end
            disp(['Writing files for ' dy '-' datestr(TS,'mmm-yyyy')]);
            
            % Make output directory if it exists
            if ~exist([outdir filesep name filesep yr filesep mo filesep dy],'file')
                mkdir([outdir filesep name filesep yr filesep mo filesep dy]);
            end
            
            % Write each variable to output file for every hour
            for hour = 0:23 
                c = c+1;
                hr = num2str(hour);
                if numel(hr) == 1, hr = ['0' hr]; end    
                VarList = {'PRES','TMP','UGRD','VGRD','SPFH','APCP','DSWRF','DLWRF'};
                VarAttList = {'Pa','K','m/s','m/s','kg/kg','kg/m^2','W/m^2','W/m^2'};
                ncfile = [outdir filesep name filesep yr filesep mo filesep dy filesep name '.' yr mo dy hr '00.nc'];
                if exist(ncfile,'file')
                    delete(ncfile);
                end
                for i = 1:numel(VarList)
                    switch i
                        case 1, var = pres(:,:,c);
                        case 2, var = airt_corr(:,:,c)+273.15;
                        case 3, var = uwind(:,:,c);
                        case 4, var = vwind(:,:,c);
                        case 5, var = spfh(:,:,c);
                        case 6, var = apcp_corr(:,:,c);
                        case 7, var = dswrf(:,:,c);
                        case 8, var = dlwrf(:,:,c);
                    end
                    nccreate(ncfile,VarList{i},'Dimensions', {'lon', sz(1),'lat', sz(2)},'DataType','single','DeflateLevel',9);
                    ncwriteatt(ncfile,['/' VarList{i}],'units',VarAttList{i});
                    ncwrite(ncfile,VarList{i},var);
                end
                
                nccreate(ncfile,'lat','Dimensions', {'lat', sz(2)},'DataType','single');
                ncwriteatt(ncfile,'/lat','units','degrees_north');
                ncwrite(ncfile,'lat',lat(locy));

                nccreate(ncfile,'lon','Dimensions', {'lon', sz(1)},'DataType','single');
                ncwriteatt(ncfile,'/lon','units','degrees_east');
                ncwrite(ncfile,'lon',lon(locx));
            end
        end
    end
end

% Remove temporary directory
if exist(tmp_dir,'file')
    rmdir(tmp_dir,'s');
end
