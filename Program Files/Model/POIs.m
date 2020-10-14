function [SpatialData] = POIs(Info, SpatialData)
% Get domains for points of interest (likely to be snow measurement 
% locations, but could potentially be other spatial features, such as 
% watersheds)
%
% Inputs: Info - Structure containing information about the model run
%         SpatialData - Structure containing spatial data for the model run
% Outputs: SpatialData - Updated structure containing spatial data
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated 04/2017.

disp(['Getting Points of Interest for ' Info.NameIdentifier]);
a = dir(Info.POIDir);
c = 0;
% Get Distributed Forcing Parameters (if any)
forcingCorrectionData = [];
if exist([Info.POIDir '/POIForcingCorrectionFactors.csv'],'file')
    fid = fopen([Info.POIDir '/POIForcingCorrectionFactors.csv']);
    tline = fgetl(fid);
    headers = strsplit(tline,',');
    tline = fgetl(fid);
    while ischar(tline)       
        forcingCorrectionData = [forcingCorrectionData; strsplit(tline,',')];
        tline = fgetl(fid);
    end
    fclose(fid);
    SpatialData.POIForcingCorrectionFactors = 1;
else
    SpatialData.POIForcingCorrectionFactors = 0;
end

% Loop through the directory containing the points of interest shapefiles
for i = 1:numel(a)
    nc = numel(a(i).name);
    if numel(a(i).name) > 4
        if strcmp(a(i).name(nc-3:nc),'.shp')
            fname = [tempname '.tif'];
            tmpshpdir = tempname;
            DEMName = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'DEM.tif'];
            mkdir(tmpshpdir);
            % Make sure that the shapefile is in the target model projection
            evalc(['!ogr2ogr -t_srs "' Info.Proj4String '" -dim 2 "' fullfile(tmpshpdir,a(i).name) '" "' fullfile(Info.POIDir,a(i).name) '"']);
            % Read the shapefile, loop through the features (which must have a "Name" attribute) <- INPORTANT
            % Instead of shaperead, use a simple function which provides the
            % relavent info without the mapping toolbox
            % S = shaperead(fullfile(tmpshpdir,a(i).name));
            [geo,S] = shapeinfo_sp([Info.ProgramFilesDir '/python'],fullfile(tmpshpdir,a(i).name));
            for j = 1:numel(S)
                c = c+1;
                imwrite2tif(zeros(size(SpatialData.z)),[],fname,'logical');
                evalc(['!python "' Info.ProgramFilesDir filesep 'Python' filesep 'gdalcopyproj.py" "' DEMName '" "' fname '"']);
                % Get a mask for each feature
                evalc(['!gdal_rasterize -at -sql "SELECT * from ' strrep(a(i).name,'.shp','') ' where Name=''' S(j).Name '''" -b 1 -burn 1 "' fullfile(tmpshpdir,a(i).name) '" "' fname '"']);
                wsmask = imread(fname);
                delete(fname)
                SpatialData.DispLocs(c).pos = find(wsmask>0);
                SpatialData.DispLocs(c).weights = ones(size(SpatialData.DispLocs(c).pos));
                SpatialData.DispLocs(c).Name = S(j).Name;
                
                SpatialData.DispLocs(c).PrecipMult = [];
                SpatialData.DispLocs(c).TSnow = [];
                SpatialData.DispLocs(c).TempAdjustment = [];
                SpatialData.DispLocs(c).X = (geo.BoundingBox(1) + geo.BoundingBox(2))/2;
                SpatialData.DispLocs(c).Y = (geo.BoundingBox(3) + geo.BoundingBox(4))/2;
                if SpatialData.POIForcingCorrectionFactors
                    for k = 1:numel(headers)
                        if strcmp(headers{k},S(j).Name)
                            SpatialData.DispLocs(c).PrecipMult = str2double(forcingCorrectionData{1,k});
                            SpatialData.DispLocs(c).TSnow = str2double(forcingCorrectionData{2,k});
                            SpatialData.DispLocs(c).TempAdjustment = str2double(forcingCorrectionData{3,k});
                        end
                    end
                end
            end
            delete([tmpshpdir filesep '*'])
            rmdir(tmpshpdir);
        end
    end
end
