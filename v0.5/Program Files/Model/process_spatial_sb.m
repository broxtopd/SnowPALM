function [SpatialData] = process_spatial_sb(Info,SpatialData,modelpars)
% Function to compute Sb Index for individual tiles
% 
% Inputs: Info is the information structure that is used throughout the snow model
%         SpatialData - Structure containing spatial data used in the model
%         modelpars - Structure containing model parameters
% Outputs: SpatialData - Updated structure containing spatial data 
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated April 2017

% Angle divisions in alphamat
Alphas = linspace(0,2*pi,Info.NAngleDivisions+1);

% Find the minimum and maximum direction (mindir, maxdir) (look 15 degrees 
% to the right and left of the specified wind direction)
PrevailingWindAngle = modelpars.PrevailingWindAngle+90;
if PrevailingWindAngle > 360
    PrevailingWindAngle = PrevailingWindAngle - 360;
end
PrevailingWindAngle_pi12 = (PrevailingWindAngle) * (pi/180) + pi/12;
PrevailingWindAngle_m_pi12 = (PrevailingWindAngle) * (pi/180)- pi/12;
if PrevailingWindAngle_pi12 > 2*pi
    PrevailingWindAngle_pi12 = PrevailingWindAngle_pi12 - 2*pi;
end
if PrevailingWindAngle_m_pi12 < 0
    PrevailingWindAngle_m_pi12 = PrevailingWindAngle_m_pi12 + 2*pi;
end
if PrevailingWindAngle_pi12 < PrevailingWindAngle_m_pi12
    mindir = find(Alphas > PrevailingWindAngle_m_pi12,1,'first');
    maxdir = find(Alphas < PrevailingWindAngle_pi12,1,'first');
else
    mindir = find(Alphas > PrevailingWindAngle_m_pi12,1,'first');
    maxdir = find(Alphas < PrevailingWindAngle_pi12,1,'last');
end

sxfile = [Info.SpatialDataDir filesep Info.SxSource '_Sx' filesep 'Sx_' num2str(mindir) '_' num2str(maxdir) '_' num2str(modelpars.Sx_Exponent) '_' num2str(modelpars.Sx_Rescale) '.tif'];
sbfile = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'Sb_' num2str(mindir) '_' num2str(maxdir) '_' num2str(modelpars.Sx_Exponent) '_' num2str(modelpars.Sx_Rescale) '_' num2str(modelpars.Sb_Sepdist) '.mat'];
% Create if the Sx file exists and the Sb file does not exist
if exist(sxfile,'file')
    if ~exist(sbfile,'file')
        % Maximum distance to look at Sx (can't include less than 3 pixels)
        IndexDist = max(3*Info.NewResolution,modelpars.Sb_Sepdist);
        xmin = num2str(Info.Cutout.min_x0-IndexDist);
        ymax = num2str(Info.Cutout.max_y0+IndexDist);
        xmax = num2str(Info.Cutout.max_x0+IndexDist);
        ymin = num2str(Info.Cutout.min_y0-IndexDist);
        tr = num2str(Info.NewResolution);

        % Read the Sx file
        fname = [tempname '.tif'];
        evalc(['!gdalwarp -dstnodata -999 -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' sxfile '" "' fname '"']);
        sx = imread(fname);
        delete(fname);
        sx(sx < 0) = 1;

        Alphas = [SpatialData.Alphas(mindir) SpatialData.Alphas(maxdir)]+pi/2;

        % Define maximum and minimum X and Y to search
        [nr,nc] = size(sx);
        min_x = round(IndexDist/Info.NewResolution+1);
        max_x = round(nc - IndexDist/Info.NewResolution);
        min_y = round(IndexDist/Info.NewResolution+1);
        max_y = round(nr - IndexDist/Info.NewResolution);
        
        % Define mask (if we are cutting out with shapefile, then don't operate on these cells) 
        mask_i = single(zeros(size(sx)));
        mask_i(min_y:max_y,min_x:max_x) = 1;
              
        if ~isfield(Info, 'ExtentFile')
            mask = single(zeros(size(sx)));
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
        mask = mask .* mask_i;

        % We will consider individual portions of the image at a time to make
        % sure that the search is not too large
        divmap = uint16(zeros(round(sqrt(sqrt(numel(sx))))));
        for i = 1:numel(divmap)
            divmap(i) = i;
        end
        divmap = imresize_sp(divmap,size(sx),'nearest');

        cs = Info.NewResolution;
        x = single((1:numel(sx(1,:))) * cs);
        x = single(repmat(x,[numel(sx(:,1)) 1]));
        y = single((1:numel(sx(:,1)))' * cs);
        y = single(repmat(y,[1 numel(sx(1,:))]));
        y = flipud(y);

        nanlocs = isnan(sx);
        divmap = divmap(:);     divmap(nanlocs) = [];
        x = x(:);               x(nanlocs) = [];
        y = y(:);               y(nanlocs) = [];
        sxmat = sx(:);          sxmat(nanlocs) = [];
        mask = mask(:);         mask(nanlocs) = [];
        sxmat = int8(100 * sxmat);
        locs = single(1:numel(sx));
        locs(nanlocs) = [];

%     try
%         [sbmat] = sbfun_mex(x,y,mask,sxmat,Alphas,divmap,IndexDist);
%     catch
%         [sbmat] = sbfun(x,y,mask,sxmat,Alphas,divmap,IndexDist);
%     end

        [sbmat] = sbfun(x,y,mask,sxmat,Alphas,divmap,IndexDist);

        % Put back into matrix form
        sb = int8(NaN([size(sx) numel(Alphas) - 1]));
        for i = 1:numel(Alphas) - 1   
            tmp = sbmat(:,i);
            Tmp = NaN(size(sx));
            Tmp(locs) = tmp;
            sb(:,:,i) = Tmp;
        end
        sb =  max(0,imresize_sp(sb(min_y:max_y,min_x:max_x),size(SpatialData.z)));

        save(sbfile,'sb')
    else
        load(sbfile);
    end

    SpatialData.Sb = single(sb) / 100;
else
    error('Cannot calculate Sb because Sx file at resolution specified by Info.SxSource is not found');
end
