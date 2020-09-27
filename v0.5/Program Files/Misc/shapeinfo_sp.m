function [geo,S] = shapeinfo_sp(pythondir,ifname)
% Use the python gdal bindings to provide relevent georeferencing
% and field name information needed for SnowPALM (python program generates 
% a text file, which is read here)

ofname = [tempname '.txt'];
eval(['!python "' pythondir '/getshapeinfo.py" "' ifname '" "' ofname '"']);
fid = fopen(ofname);
ulx = fgetl(fid);
uly = fgetl(fid);
lrx = fgetl(fid);
lry = fgetl(fid);

S = [];
i = 0;
tline = fgetl(fid);
while ischar(tline)
    i = i+1;
    S(i).Name = tline;
    tline = fgetl(fid);
end
fclose(fid);
delete(ofname);
geo.BoundingBox(1) = str2double(ulx);
geo.BoundingBox(2) = str2double(lrx);
geo.BoundingBox(3) = str2double(lry);
geo.BoundingBox(4) = str2double(uly);
