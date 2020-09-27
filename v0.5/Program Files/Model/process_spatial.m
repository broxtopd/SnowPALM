function [SpatialData] = process_spatial(Info,modelpars)
% Function to load necessary spatial data and to data for future use in the
% model
%
% Inputs: Info - Structure containing information about the model run
%         modelpars - Structure containing model parameters
% Outputs: SpatialData - Structure containing spatial data used in the model
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated July 2017

disp(['Getting spatial data for ' Info.NameIdentifier]);

% Check to see if a state file already exists for this area.  If not, then
% process the spatial data.  If so, load the file.  If you want to delete
% the file and start over, delete it from the 'Info.ModelStatesDir/Spatial'
% directory, or use the 'remove_files.m' script
SpatialData = [];
NoDataVal = -9999;  % Anything less than this value will be considered NaN

% Get the Bounding boxes for each tile
xmin = num2str(Info.Cutout.min_x0);
ymax = num2str(Info.Cutout.max_y0);
xmax = num2str(Info.Cutout.max_x0);
ymin = num2str(Info.Cutout.min_y0);
% ...and the model resolution
tr = num2str(Info.NewResolution);

% Make the directory to store the spatial information
if ~exist([Info.SpatialDataDir filesep Info.NameIdentifier],'file')
    mkdir([Info.SpatialDataDir filesep Info.NameIdentifier])
end

nodatafile = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'NoData.txt'];
if ~exist(nodatafile,'file')
    % File in which to store the spatial data
    fname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'SpatialData.mat'];

    % If it doesn't exist, first clip the DTM for the tile (also clip using the
    % shapefile, if one is specified
    ofname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'DEM.tif'];
    ofname2 = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'DEM_mask.tif'];
    if ~exist(ofname2,'file'); 
        eval(['!gdalwarp -overwrite -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.DTMFile '" "' ofname '"']);
        if isfield(Info, 'ExtentFile')
            tmpshpdir = tempname;
            mkdir(tmpshpdir);
            evalc(['!ogr2ogr -t_srs "' Info.Proj4String '" -dim 2 -clipdst ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' fullfile(tmpshpdir,'tmp.shp') '" "' Info.ExtentFile '"']);
            S = shapeinfo_sp([Info.ProgramFilesDir '/python'],fullfile(tmpshpdir,'tmp.shp'));
            z = double(imread(ofname));
            if (S.BoundingBox(1)) < 1E30
                evalc(['!gdal_rasterize -ot Byte -init 0 -burn 1 -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' fullfile(tmpshpdir,'tmp.shp') '" "' ofname2 '"']);
            else
                imwrite2tif(zeros(size(z)),[],ofname2,'single');
                evalc(['!python "' Info.ProgramFilesDir '/python/gdalcopyproj.py" "' ofname '" "' ofname2 '"']);
            end
            rmdir(tmpshpdir,'s')
        end
    end

    z_0 = double(imread(ofname));
    % Instead of geotiffinfo, use a simple function which provides the
    % relavent info without the mapping toolbox
%     geo = geotiffinfo(ofname);
    geo = geotiffinfo_sp([Info.ProgramFilesDir '/python'],ofname);
    
    SpatialData.z = z_0;
    if exist(ofname2,'file');
        mask = imread(ofname2);
    else
        mask = ones(size(z_0));
    end

    nodatalocs = SpatialData.z <= NoDataVal | isnan(SpatialData.z);
    nanlocs = SpatialData.z <= NoDataVal | isnan(SpatialData.z) | mask == 0;

    if sum(~nanlocs(:)) == 0
        fid = fopen(nodatafile,'w');
        fwrite(fid,'Nothing to be done');
        fclose(fid);
        delete(ofname);
        if exist(ofname2,'file')
            delete(ofname2);
        end
        SpatialData = [];
    end
end

ofname2 = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'DEM_mask.tif'];
if exist(nodatafile,'file') && exist(ofname2,'file')
    delete(ofname2);
end
    
    
if ~exist(nodatafile,'file')
    if ~exist(fname,'file') 

        %% For Tilting Experiments
    % %     % Flat
    % %     SpatialData.z(:) = mean(z_0(:));
    %     % South Facing
    %     SpatialData.z = 0.3 * (repmat(linspace(str2double(ymax),str2double(ymin),numel(z_0(:,1)))',[1 numel(z_0(1,:))]) - str2double(ymin)) + mean(z_0(:));
    % %     % North Facing
    % %     SpatialData.z = 0.3 * (repmat(linspace(str2double(ymin),str2double(ymax),numel(z_0(:,1)))',[1 numel(z_0(1,:))]) - str2double(ymin)) + mean(z_0(:));
    %     % pcolor(flipud(double(SpatialData.z)))
    %     % shading flat
    %     % colorbar
    %     % pause
        %%

        % Get the DEM for the forcing interpolation
        SpatialData.z(nodatalocs) = mean(SpatialData.z(~nodatalocs));

        % Matrices for x and y coordinants (local coordinate system)
        SpatialData.cellsize = sqrt(abs(geo.RefMatrix(2)) * abs(geo.RefMatrix(4)));
        SpatialData.x = repmat((geo.BoundingBox(1)+SpatialData.cellsize/2:SpatialData.cellsize:geo.BoundingBox(2)-SpatialData.cellsize/2),[numel(SpatialData.z(:,1)),1]);
        SpatialData.y = repmat((geo.BoundingBox(4)-SpatialData.cellsize/2:-SpatialData.cellsize:geo.BoundingBox(3)+SpatialData.cellsize/2)',[1,numel(SpatialData.z(1,:))]);

        % Find the slope and aspect
        [fx,fy] = gradient(SpatialData.z,SpatialData.cellsize,SpatialData.cellsize); % uses simple, unweighted gradient of immediate neighbours
        [R,S]=cart2pol(fy,fx);  % convert to cartesian coordinates
        S=single(atan(S));    	
        R=single(R.*-1+pi);     % convert asp 0 facing south
        SpatialData.S = S;
        R2 = R + 2*pi; R(R<0) = 0; R2(R>0) = 0; R = R+R2;
        SpatialData.R = R;  
        nanlocs = nanlocs | SpatialData.R <= NoDataVal;

        % Get matrices of longitude and latitude coordinates (for forcing interpolation)
        ifname = ofname;
        ofname1 = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'lat.tif'];
        ofname2 = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'lon.tif'];
        if exist(ofname1,'file')
            delete(ofname1);
        end
        if exist(ofname2,'file');
            delete(ofname2)
        end
        evalc(['!python "' Info.ProgramFilesDir filesep 'Python' filesep 'pixels2latlonmat.py" "' ifname '" "' ofname1 '" "' ofname2 '"']);
        SpatialData.LatMap = double(imread(ofname1));
        SpatialData.LonMap = double(imread(ofname2));
        delete(ofname1);
        delete(ofname2);

        % Get the vegetation height map
        ofname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'Veght.tif'];
        if ~exist(ofname,'file'); 
            evalc(['!gdalwarp -overwrite -r average -t_srs "' Info.Proj4String '" -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.VegHeightFile '" "' ofname '"']);
        end
        SpatialData.canopy_height = double(imread(ofname));
        delete(ofname);
        SpatialData.canopy_height(SpatialData.canopy_height < 0) = 0;
        SpatialData.canopy_height(isnan(SpatialData.canopy_height)) = 0;
        nanlocs = nanlocs | SpatialData.canopy_height <= NoDataVal;

        % Get the vegetation cover map
        ofname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'Cover.tif'];
        if ~exist(ofname,'file'); 
            evalc(['!gdalwarp -overwrite -r average -t_srs "' Info.Proj4String '" -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.VegCoverFile '" "' ofname '"']);
        end
        SpatialData.canopy_coverage = double(imread(ofname));
        delete(ofname);
        SpatialData.canopy_coverage = max(0,min(100,SpatialData.canopy_coverage));
        nanlocs = nanlocs | SpatialData.canopy_height <= NoDataVal;

        SpatialData.z = z_0;

        %% Tilting Experiments
    %     ymax = num2str(Info.NSWE(1));
    %     ymin = num2str(Info.NSWE(2));
    %     xmin = num2str(Info.NSWE(3));
    %     xmax = num2str(Info.NSWE(4));
    %     ofname = [tempname '.tif'];
    %     evalc(['!gdalwarp -overwrite -t_srs "' Info.Proj4String '" -r bilinear -tr 100 100 -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.DTMFile '" "' ofname '"']);
    %     im = imread(ofname);
    %     SpatialData.z(:) = mean(im(:));
    %     SpatialData.z_forcing(:) = mean(im(:));
    %     delete(ofname);
        %%

        SpatialData.z(nanlocs) = NaN;
        SpatialData.S(nanlocs) = NaN;
        SpatialData.R(nanlocs) = NaN;
        SpatialData.canopy_coverage(nanlocs) = NaN;
        SpatialData.canopy_height(nanlocs) = NaN;

        save(fname,'SpatialData');
    else
        load(fname);
    end

    if Info.ComputeAngles
        % Find horizon angles (without vegetation)
        SpatialData = process_spatial_angles(Info,SpatialData,0); 
        % Find closeness to vegetation and compute closeness indices
        [SpatialData] = process_spatial_candist(Info,SpatialData);
        rad_idx_shortwave = round(modelpars.EFLength) == SpatialData.dstpars;
        SpatialData.WgtFun_shortwave = SpatialData.WgtFun(:,:,rad_idx_shortwave);
        SpatialData.S_0 = SpatialData.S;
        SpatialData.R_0 = SpatialData.R;
        z = SpatialData.z;
        z(isnan(z)) = mean(z(~isnan(z)));
        [fx,fy] = gradient(SpatialData.WgtFun_shortwave+z,SpatialData.cellsize,SpatialData.cellsize); % uses simple, unweighted gradient of immediate neighbours
        [R,S]=cart2pol(fy,fx); % convert to cartesian coordinates
        S=single(atan(S));    % steepest slope
        R=single(R.*-1+pi);     % convert asp 0 facing south
        R(isnan(SpatialData.z)) = NaN;
        S(isnan(SpatialData.z)) = NaN;
        SpatialData.S_veg = S;
        R2 = R + 2*pi; R(R<0) = 0; R2(R>0) = 0; R = R+R2;
        SpatialData.R_veg = R;
        % Find horizon angles (with vegetation)
        SpatialData = process_spatial_angles(Info,SpatialData,1); 
    end

    % Find wind indices
    
    if Info.ComputeSxIndex
        SpatialData = process_spatial_sx(Info,SpatialData,modelpars);
    else
        SpatialData.Sx = ones(size(SpatialData.z));
        SpatialData.Sx(isnan(SpatialData.z)) = NaN;
    end
    if Info.ComputeSbIndex
        SpatialData = process_spatial_sb(Info,SpatialData,modelpars); 
    else
        SpatialData.Sb = zeros(size(SpatialData.z));
        SpatialData.Sb(isnan(SpatialData.z)) = NaN;
    end
end
