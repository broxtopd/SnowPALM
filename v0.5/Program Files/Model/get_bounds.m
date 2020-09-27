function [Info,Cutout] = get_bounds(Info)
% Returns the bounding box of each model tile

% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated 07/2017
%%
% MaxSquarePixles contains the maximum number of pixels in a tile, so the
% side of a tile will have up to the square root of that
boxside = round(sqrt(Info.MaxSquarePixels));

% If the bounds of the model subdomin are determined using a vector data file... 
if ~isfield(Info,'NSWE')
    tmpshpdir = tempname;
    mkdir(tmpshpdir);
    % Get the bounding box for the shapefile
    evalc(['!ogr2ogr -t_srs "' Info.Proj4String '" -dim 2 "' fullfile(tmpshpdir,'tmp.shp') '" "' Info.ExtentFile '"']);
    % Instead of shapeinfo, use a simple function which provides the
    % relavent info without the mapping toolbox
    % shpinfo = shapeinfo(fullfile(tmpshpdir,'tmp.shp'));
  	[shpinfo,~] = shapeinfo_sp([Info.ProgramFilesDir '/python'],fullfile(tmpshpdir,'tmp.shp'));
    
    min_x = shpinfo.BoundingBox(1)-1;
    max_x = shpinfo.BoundingBox(2)+1;
    min_y = shpinfo.BoundingBox(3)-1;
    max_y = shpinfo.BoundingBox(4)+1;
    
    delete([tmpshpdir filesep '*'])
 	rmdir(tmpshpdir);
    Info.NSWE = [max_y min_y min_x max_x];
% Otherwise, if they are prescribed
else
    min_x = round(Info.NSWE(3));
    max_x = round(Info.NSWE(4));
    min_y = round(Info.NSWE(2));
    max_y = round(Info.NSWE(1));
end

% Set up vectors defining the discretization of each model tile...
c = min_x:boxside*Info.NewResolution:max_x-1;
c = [c max_x];
r = min_y:boxside*Info.NewResolution:max_y-1;
r = [r max_y];

% ...and save the model bounds of each model tile to an output structure
k = 0;
for i = 1:numel(c) - 1
    for j = 1:numel(r) - 1
        k = k+1;
        Cutout(k).min_x0 = c(i);
        Cutout(k).max_x0 = c(i+1);
        Cutout(k).min_y0 = r(j);
        Cutout(k).max_y0 = r(j+1);     
    end
end
