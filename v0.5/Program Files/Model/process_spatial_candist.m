function [SpatialData] = process_spatial_candist(Info,SpatialData)
% Function to process the angles-to-horizon used to find shadows in the
% solar forcing program.  
% 
% Inputs: Info is the information structure that is used throughout the snow model
%         SpatialData - Structure containing spatial data used in the model
% Outputs: SpatialData - Updated structure containing spatial data 
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated November 2012

NoDataVal = -9999;      % Anything less than this value will be considered NaN
dstpars = [1 2 3 4 5 6 7 8 9 10];         % Areas of influence (1-10m)

fname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'CanIdx.mat'];
% Check to see if a state file already exists for this area.  If not, then
% process the spatial data.  If so, load the file.  
if ~exist(fname,'file'); 
    
    maxdst = Info.CanDist / Info.NewResolution; 	% Maximum distance to search for obstructions [M]

    xmin = num2str(Info.Cutout.min_x0);
    ymax = num2str(Info.Cutout.max_y0);
    xmax = num2str(Info.Cutout.max_x0);
    ymin = num2str(Info.Cutout.min_y0);
    tr = num2str(Info.NewResolution);

    % Load the vegetation cover fi
    evalc(['!gdalwarp -r average -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.VegCoverFile '" "tmp_' Info.NameIdentifier '.tif"']);
    vegcover = double(imread(['tmp_' Info.NameIdentifier '.tif']));
    delete(['tmp_' Info.NameIdentifier '.tif']);
    vegcover = max(0,100-vegcover);

    evalc(['!gdalwarp -r average -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' Info.VegHeightFile '" "tmp_' Info.NameIdentifier '.tif"']);
    veght = double(imread(['tmp_' Info.NameIdentifier '.tif']));
    delete(['tmp_' Info.NameIdentifier '.tif']);
    cs = Info.NewResolution; 
    [nr,nc] = size(veght);

    % Find the vegetation distances (put in separate subfunction which can
    % be compiled because it is quite slow)
    disp('Finding distances to Vegetation:');

    % Set up x an y matricies
    x = single((1:numel(vegcover(1,:))) * cs);
    x = single(repmat(x,[numel(vegcover(:,1)) 1]));
    y = single((1:numel(vegcover(:,1)))' * cs);
    y = single(repmat(y,[1 numel(vegcover(1,:))]));
    y = flipud(y);

    % We will consider individual portions of the image at a time to make
    % sure that the search is not too large
    divmap = uint16(zeros(round(sqrt(sqrt(numel(vegcover))))));
    for i = 1:numel(divmap)
        divmap(i) = i;
    end

    % Here, we define smaller chunks to operate on
    divmap = imresize_sp(divmap,size(vegcover),'nearest');
    
    % We need a different mask because the algorithm should look beyond the
    % cutout mask by 'maxdist' meters
    
    if ~isfield(Info, 'ExtentFile')
        mask = single(ones(size(vegcover)));
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
    
    nanlocs = vegcover <= NoDataVal;
    divmap = divmap(:);     divmap(nanlocs) = [];
    x = x(:);               x(nanlocs) = [];
    y = y(:);               y(nanlocs) = [];
    vegcover = vegcover(:);	vegcover(nanlocs) = [];
    veght = veght(:);     	veght(nanlocs) = [];
    mask = mask(:);         mask(nanlocs) = [];
    locs = single(1:numel(vegcover));


    % Find the distances (put in separate subfunction which can be 
    % converted to c and compiled because it is quite slow)

    tic
%     try
%         [canpos,wgtfun] = dstfun_mex(x,y,vegcover,veght,mask,divmap,maxdst,dstpars);
%     catch
        [canpos,wgtfun] = dstfun(x,y,vegcover,veght,mask,divmap,maxdst,dstpars);
%     end
    toc 


    % Put back into matrix form
    CanPos = single(NaN([nr,nc]));
    for j = 1:numel(canpos)
        CanPos(locs(j)) = canpos(j);
    end
    WgtFun = zeros([[nr,nc] numel(dstpars)]);
    for i = 1:numel(dstpars)
        tmp = single(NaN([nr,nc]));
        for j = 1:numel(canpos)
            tmp(locs(j)) = wgtfun(i,j);
        end
        WgtFun(:,:,i) = tmp;
    end

    % Save the results to a file
    save(fname,'CanPos','WgtFun','dstpars');
else
    load(fname);
end

SpatialData.CanPos = CanPos;
SpatialData.WgtFun = WgtFun;
SpatialData.dstpars = dstpars;

