function combine_sx_maps(Info,IdentifierList,modelpars)
% Combines all of the map tiles from the model chunk data created by the display_output program
%
% Inputs: Info - Structure containing information about the model run
%         IdentifierList - List of map tiles to combine
%         modelpars - Structure containing values of the model parameters
% No outputs, instead generates GIS files
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

% Folder and name given to the output file
fpath = [Info.SpatialDataDir filesep Info.NameIdentifier '_Sx'];
fname = ['Sx_' num2str(mindir) '_' num2str(maxdir) '_' num2str(modelpars.Sx_Exponent) '_' num2str(modelpars.Sx_Rescale)];

% If this file does not exist ...
if ~exist([fpath filesep fname '.tif'],'file')
    display('Combining Sx Maps...')
    % Combine the maps, create an output tif
    if ~exist([fpath filesep fname '.tif'],'file')
        combine_vrt_tif(IdentifierList,fname,fpath);
        for i = 1:numel(IdentifierList)
            try
                rmdir([fpath filesep IdentifierList(i).Identifier],'s');
            end
        end
    end
end

function combine_vrt_tif(IdentifierList,fname,fpath)
% Builds vrt's out of all of the tiles generated for each model tile, and
% then converts to output GIS format (currently geotiff, but could be
% modified to be any / all gdal supported formats)

% Create a list of all of the files to put together
flist = [tempname '.txt'];
fid = fopen(flist,'w');
for i = 1:numel(IdentifierList)
    addr = [fpath filesep IdentifierList(i).Identifier filesep fname '.tif'];
    if exist(addr,'file')
        fprintf(fid,[strrep(addr,'\','/') '\n']);
    end
end
fclose(fid);

% Create the folder where we will put the data if it does not exist
if ~exist(fpath,'file')
    mkdir(fpath);
end

% Create the tiled mosaic and convert to geotiff
eval(['!gdalbuildvrt -input_file_list "' flist '" "' fpath filesep fname '.vrt"']);
eval(['!gdal_translate "' fpath filesep fname '.vrt" "' fpath filesep fname '.tif"']);
% The following code adds overview tiles (for fast GIS display of large datasets) - disabled for now
% evalc(['!gdaladdo "' fpath filesep fname '.tif" 4 16 64']);
% Delete all of the temporary files
delete([fpath filesep fname '.vrt']);
delete(flist);
   