function [SpatialData] = process_spatial_sx(Info,SpatialData,modelpars)
% Function to compute Sx Index for individual tiles
% 
% Inputs: Info is the information structure that is used throughout the snow model
%         SpatialData - Structure containing spatial data used in the model
%         modelpars - Structure containing model parameters
% Outputs: SpatialData - Updated structure containing spatial data 
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated April 2017

cover = SpatialData.canopy_coverage / 100;
slope = sin(SpatialData.S);
Alphas = linspace(0,2*pi,Info.NAngleDivisions+1);   % Angle divisions in alphamat

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

% Get list of directions to look for
if maxdir >= mindir
    dirs = mindir:maxdir;
else
    dirs = [mindir:numel(Alphas)-1 1:maxdir];
end

% Output map names
fpath = [Info.SpatialDataDir filesep Info.SxSource '_Sx'];
% Map created here
fname = ['Sx_' num2str(mindir) '_' num2str(maxdir) '_' num2str(modelpars.Sx_Exponent) '_' num2str(modelpars.Sx_Rescale)];

SpatialData.Sx = [];
% If this file does not exist ...
if ~exist([fpath filesep fname '.tif'],'file');
    % Create a Sx Map from the horizon angles if we are operating on the
    % SxSource domain
    if strcmp(Info.SxSource,Info.NameIdentifier_orig)
        sx = (1-(mean(1-sin(double(SpatialData.VegAlphaMat(:,:,dirs))/(255/(pi/2))),3) .* (1-cover)) .^ modelpars.Sx_Exponent);
        sx_min = 1-slope;
        sx = max(sx,sx_min);
        sx = max(0,(sx - modelpars.Sx_Rescale)*1/(1-modelpars.Sx_Rescale));
        fpath = [fpath  filesep Info.NameIdentifier];
        ofname = [fpath filesep fname '.tif'];
        if ~exist(fpath,'file')
            mkdir(fpath)
        end
        % Copy georeference information from the DEM
        ifname = [Info.SpatialDataDir filesep Info.NameIdentifier filesep 'DEM.tif'];
        imwrite2tif(sx,[],ofname,'single');
        eval(['!python "' Info.ProgramFilesDir filesep 'Python' filesep 'gdalcopyproj.py" "' ifname '" "' ofname '"']);
    else
        error('Sx map at the resolution specified by Info.SxSource is not found');
    end
% else, cut the Sx map from the Sx source map
else
    xmin = num2str(Info.Cutout.min_x0);
    ymax = num2str(Info.Cutout.max_y0);
    xmax = num2str(Info.Cutout.max_x0);
    ymin = num2str(Info.Cutout.min_y0);
    tr = num2str(Info.NewResolution);
    ofname = [tempname '.tif'];
    evalc(['!gdalwarp -overwrite -t_srs "' Info.Proj4String '" -r bilinear -tr ' tr ' ' tr ' -te ' xmin ' ' ymin ' ' xmax ' ' ymax ' "' fpath filesep fname '.tif" "' ofname '"']);
    SpatialData.Sx = imread(ofname);
    delete(ofname);
end
