function [SpatialData] = process_spatial_angles(Info,SpatialData,useveg)
% Function to process the angles-to-horizon used to find shadows in the
% solar forcing program.  
% 
% Inputs: Info is the information structure that is used throughout the snow model
%         SpatialData - Structure containing spatial data used in the model
%         useveg - flag indicating whether to load the vegetation data when
%           computing angles (difference between hard shading and soft
%           shading)
% Outputs: SpatialData - Updated structure containing spatial data 
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated September 2016
%%
NoDataVal = -9999;      % Anything less than this value will be considered NaN
NAngleDivisions = Info.NAngleDivisions;          % Number of anglular divisions for the analysis (in this case, divide into 1 degree increments)
orig_sz = size(SpatialData.z);
SpatialData_angles.AlphaMat = uint8(zeros([orig_sz,NAngleDivisions]));
BaseDist = 50;   

if useveg
    fname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'VegAngles_' num2str(NAngleDivisions) '.mat'];
else
    fname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'Angles_' num2str(NAngleDivisions) '.mat'];
end

% This code runs iteratively, computing horizon angles for larger and
% larger areas at successively lower resolutions (to take into account
% distant shadowing effects)
nlevels = ceil(log(Info.MaxTerrainDist/(BaseDist*Info.NewResolution)/log(2))) + 1;

% Check to see if a state file already exists for this area.  If not, then
% process the spatial data.  If so, load the file.  
if ~exist(fname,'file'); 
    for level = 1:nlevels
        xmin = num2str(Info.Cutout.min_x0-Info.NewResolution*BaseDist*2^(level-1));
        ymax = num2str(Info.Cutout.max_y0+Info.NewResolution*BaseDist*2^(level-1));
        xmax = num2str(Info.Cutout.max_x0+Info.NewResolution*BaseDist*2^(level-1));
        ymin = num2str(Info.Cutout.min_y0-Info.NewResolution*BaseDist*2^(level-1));
        tr = num2str(Info.NewResolution*2^(level-1));

        evalc(['!gdalwarp -overwrite -r bilinear -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.DTMFile '" "tmp_' Info.NameIdentifier '.tif"']);
        dem = double(imread(['tmp_' Info.NameIdentifier '.tif']));
        delete(['tmp_' Info.NameIdentifier '.tif']);
        
        %% For Tilting Experiments
%         % Flat
%         dem(:) = mean(dem(dem > -999));
%         % South Facing
%         dem = 0.3 * (repmat(linspace(str2double(ymax),str2double(ymin),numel(dem(:,1)))',[1 numel(dem(1,:))]) - str2double(ymin)) + mean(dem(dem > -999));
%         % North Facing
%         dem = 0.3 * (repmat(linspace(str2double(ymin),str2double(ymax),numel(dem(:,1)))',[1 numel(dem(1,:))]) - str2double(ymin)) + mean(dem(dem > -999));

%         pcolor(flipud(double(dem)))
%         shading flat
%         colorbar
%         pause
        %%
  
        % If considering vegetation, load the canopy height map
        if useveg
            evalc(['!gdalwarp -r average -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.VegHeightFile '" "tmp_' Info.NameIdentifier '.tif"']);
            veght = double(imread(['tmp_' Info.NameIdentifier '.tif']));
            delete(['tmp_' Info.NameIdentifier '.tif']);
            dem = dem + veght;
        end
        [nr,nc] = size(dem);

        % mindist is the minimum distance to look for obstructions at a
        % given telescoping level
        % maxdist is the maximum distance to look for obstructions at a
        % given telescoping level
        if level == 1
            mindist = 0;
            maxdst = BaseDist*Info.NewResolution;
        else
            mindist = BaseDist*2^(level-2)*Info.NewResolution;
            maxdst = BaseDist*2^(level-1)*Info.NewResolution;
        end
        min_x = BaseDist + 1;
        max_x = nc - BaseDist;
        min_y = BaseDist + 1;
        max_y = nr - BaseDist;
        cs = Info.NewResolution*2^(level-1);

        % Initialize the Angles
        Alphas = linspace(0,2*pi,NAngleDivisions+1);
        Alphas = Alphas+pi/2;

        demlocs = single(1:numel(dem));

        % Set up x an y matricies
        x = single((1:numel(dem(1,:))) * cs);
        x = single(repmat(x,[numel(dem(:,1)) 1]));
        y = single((1:numel(dem(:,1)))' * cs);
        y = single(repmat(y,[1 numel(dem(1,:))]));
        y = flipud(y);
        z = single(dem);
        z(z<=NoDataVal) = NaN;
        

        % We will consider individual portions of the image at a time to make
        % sure that the search is not too large
        divmap = uint16(zeros(round(sqrt(sqrt(numel(z))))));
        for i = 1:numel(divmap)
            divmap(i) = i;
        end

        divmap = imresize_sp(divmap,size(z),'nearest');

        % We need a different mask because the algorithm should look beyond the
        % cutout mask by 'maxdist' meters

        if ~isfield(Info, 'ExtentFile')
            mask = single(zeros(size(z)));
            mask(min_y:max_y,min_x:max_x) = 1;
        else
            tmpshpdir = tempname;
            mkdir(tmpshpdir);
            evalc(['!ogr2ogr -t_srs "' Info.Proj4String '" -dim 2 -clipdst ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' fullfile(tmpshpdir,'tmp.shp') '" "' Info.ExtentFile '"']);
            evalc(['!gdal_rasterize -ot Byte -init 0 -burn 1 -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' fullfile(tmpshpdir,'tmp.shp') '" "tmp_' Info.NameIdentifier '_mask.tif"']);
            mask = double(imread(['tmp_' Info.NameIdentifier '_mask.tif']));
            delete(['tmp_' Info.NameIdentifier '_mask.tif']);
            rmdir(tmpshpdir,'s')
        end

        % Vectorize everything, if outside of the map area (z is NaN), do not
        % consider
        nanlocs = isnan(z);
        divmap = divmap(:);     divmap(nanlocs) = [];
        x = x(:);               x(nanlocs) = [];
        y = y(:);               y(nanlocs) = [];
        z = z(:);               z(nanlocs) = [];
        mask = mask(:);         mask(nanlocs) = [];
        demlocs = demlocs(:);   demlocs(nanlocs)   = [];

        % Perform the angle extraction (put in separate subfunction which can
        % be converted to c and compiled because it is quite slow)
        disp(['Starting Angle Extraction pass ' num2str(level)]);
        tic
        % If using compiled code (explained in function)
%         try
%             [alphamat] = anglefun_mex(x,y,z,mask,Alphas,divmap,mindist,maxdst);
%         catch
            [alphamat] = anglefun(x,y,z,mask,Alphas,divmap,mindist,maxdst);
%         end
        toc

        % Put back into matrix form
        AlphaMat = uint8(NaN([size(dem) numel(Alphas) - 1]));

        for i = 1:numel(Alphas) - 1
            tmp = alphamat(:,i);
            Tmp = NaN(size(dem));
            Tmp(demlocs) = tmp;
            
            if sum(Tmp(:)==0) == numel(Tmp(~isnan(Tmp)))
                if i > 1
                    AlphaMat(:,:,i) = AlphaMat(:,:,i-1);
                end
            else
                AlphaMat(:,:,i) = Tmp;
            end
        end
        Tmp = AlphaMat(:,:,1);
        if sum(Tmp(:)==0) == numel(Tmp(~isnan(Tmp)))
             AlphaMat(:,:,1) = AlphaMat(:,:,2);
        end

        AlphaMat = AlphaMat(min_y:max_y,min_x:max_x,:);
        for i = 1:numel(Alphas) - 1
            if ~isempty(AlphaMat)
                SpatialData_angles.AlphaMat(:,:,i) = max(SpatialData_angles.AlphaMat(:,:,i),uint8(imresize_sp(double(AlphaMat(:,:,i)),orig_sz,'bilinear')));
            end
        end
        SpatialData_angles.Alphas = single(Alphas-pi/2);
        
        clear dem;
        clear x;
        clear y;
        clear z;
        clear demlocs;
        clear alphamat;
    
    end
    save(fname,'SpatialData_angles');
else
    load(fname);
end

% Put the data into the appropriate variable, based on whether it is
% horizon angles based on ground only, or ground + veg
if useveg
    SpatialData.VegAlphaMat = SpatialData_angles.AlphaMat;
    SpatialData.VegAlphas = SpatialData_angles.Alphas;
else
    SpatialData.AlphaMat = SpatialData_angles.AlphaMat;
    SpatialData.Alphas = SpatialData_angles.Alphas;
end           