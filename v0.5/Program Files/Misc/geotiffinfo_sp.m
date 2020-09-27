function [geo] = geotiffinfo_sp(pythondir,ifname)
% Use the python gdal bindings to provide relevent georeferencing
% information needed for SnowPALM (python program generates a text file,
% which is read here)

ofname = [tempname '.txt'];
eval(['!python "' pythondir '/getgeorefinfo.py" "' ifname '" "' ofname '"']);
fid = fopen(ofname);
pixelWidth = fgetl(fid);
pixelHeight = fgetl(fid);
ulx = fgetl(fid);
uly = fgetl(fid);
lrx = fgetl(fid);
lry = fgetl(fid);
fclose(fid);
delete(ofname);
geo.RefMatrix(2) = str2double(pixelWidth);
geo.RefMatrix(4) = str2double(pixelHeight);
geo.BoundingBox(1) = str2double(ulx);
geo.BoundingBox(2) = str2double(lrx);
geo.BoundingBox(3) = str2double(lry);
geo.BoundingBox(4) = str2double(uly);
