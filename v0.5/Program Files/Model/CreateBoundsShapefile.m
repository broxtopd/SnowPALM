function CreateBoundsShapefile(Info,Cutout)
% Creates a shapefile showing where the map tiles are located
% This function requires the mapping toolbox
for i=1:numel(Cutout)
    NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
    S(i).Geometry = 'Polygon';
    S(i).Name = NameIdentifier;
    S(i).X = [Cutout(i).min_x0 Cutout(i).max_x0 Cutout(i).max_x0 Cutout(i).min_x0 Cutout(i).min_x0];
    S(i).Y = [Cutout(i).max_y0 Cutout(i).max_y0 Cutout(i).min_y0 Cutout(i).min_y0 Cutout(i).max_y0];
end

fpath = [Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier];
if ~exist(fpath,'file')
    mkdir(fpath);
end
shapewrite(S,[fpath filesep 'Index_tmp.shp'])
delete([fpath filesep 'Index.*']);
eval(['!ogr2ogr -a_srs "' Info.Proj4String '" "' fpath filesep 'Index.shp" "' fpath filesep 'Index_tmp.shp"']);
delete([fpath filesep 'Index_tmp.*']);

