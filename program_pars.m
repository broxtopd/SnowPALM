function [Info] = program_pars(SubDomain,Info)

%% Paths
Info.ProgramFilesDir = [pwd '/Program Files'];                                  % Required - Directory where model program files are located
Info.ModelStatesDir = [pwd '/Model States'];                                    % Required - Directory where model states are stored
Info.SpatialDataDir = [pwd '/Spatial Data'];                                    % Required - Directory where processed model spatial data are stored
Info.DisplayDir = [pwd '/Display'];                                             % Required - Directory where the model outputs (tabular data and 
                                                                                %   GIS files) are stored
Info.POIDir = [pwd '/POIs'];                                                    % Required - Directory where model points of interest (which will appear
                                                                                %   in the tabular model output) are stored
Info.ForcingDirs{1} = [pwd '/Forcing Data/prism_corrected_nldas_valles'];       % Required - Directory containing the processed forcing data files
                                                                                %   (run the get forcing program first)
% Info.ForcingDirs{2} = ['...'];    % Can have different sources for different variables (specified below)

%% Spatial Information

Info.Proj4String = 'EPSG: 26913';                                  % Model projection (does not have to match that of input 
                                                                                        %   spatial data files)
Info.DTMFile = 'GIS Files/jrb_dtm.tif';                	% Bare earth DTM file 
Info.VegHeightFile = 'GIS Files/jrb_veght.tif';          % Vegetation height file
Info.VegCoverFile = 'GIS Files/jrb_cover.tif';           % Vegetation coverage file, typically computed as canopy closure

% Model Subdomains
if strcmp(SubDomain,'Test')
    Info.NewResolution = 1;
    Info.NSWE = [3972800 3971800 361000 362000];
%     Info.ExtentFile = 'PATH TO FILE';     % File containing a polygon that contains the extent of the subdomain
    Info.SxSource = 'Test';
    Info.MaxSquarePixels = 250^2;
else
    error(['SubDomain named " ' SubDomain '" has no associated spatial information!']);
end

%% Forcing Information

Info.Forcing_ULLR = [-107 36.5 -106 35.5];                                              % Required - Bounding box that forcing data is extracted for (must be
                                                                                        %   large enough to capture enough NLDAS pixels outside of the model 
                                                                                        %   subdomains to compute lapse rates, but smaller domains are faster)

% Specify how different forcing variables are treated                                                                                         
Info.Forcing.var{1} = 'PRES';                                                           % Variable name
Info.Forcing.source{1} = 1;                                                             % Source for each forcing variable (sources are specified above)
Info.Forcing.interpMethod{1} = 'linear';                                                % Interpolation type (nearest, linear, cubic, spline, etc)
Info.Forcing.useLapseRate{1} = 1;                                                       % 1) Use a lapse rate (i.e. interpolate lapse rate variables) or 0) do simple interpolation

Info.Forcing.var{2} = 'TMP';
Info.Forcing.source{2} = 1;
Info.Forcing.interpMethod{2} = 'linear';
Info.Forcing.useLapseRate{2} = 1;

Info.Forcing.var{3} = 'UGRD';
Info.Forcing.source{3} = 1;
Info.Forcing.interpMethod{3} = 'linear';
Info.Forcing.useLapseRate{3} = 1;

Info.Forcing.var{4} = 'VGRD';
Info.Forcing.source{4} = 1;
Info.Forcing.interpMethod{4} = 'linear';
Info.Forcing.useLapseRate{4} = 1;

Info.Forcing.var{5} = 'SPFH';
Info.Forcing.source{5} = 1;
Info.Forcing.interpMethod{5} = 'linear';
Info.Forcing.useLapseRate{5} = 1;

Info.Forcing.var{6} = 'APCP';
Info.Forcing.source{6} = 1;
Info.Forcing.interpMethod{6} = 'linear';
Info.Forcing.useLapseRate{6} = 1;

Info.Forcing.var{7} = 'DLWRF';
Info.Forcing.source{7} = 1;
Info.Forcing.interpMethod{7} = 'linear';
Info.Forcing.useLapseRate{7} = 1;

Info.Forcing.var{8} = 'DSWRF';
Info.Forcing.source{8} = 1;
Info.Forcing.interpMethod{8} = 'linear';
Info.Forcing.useLapseRate{8} = 1;


%% Additional Options (move inside the subdomains if you want each subdomain simulation to have different parameters)

Info.CanDist = 50;                      % Maximum distance to consider canopy area of influence effects (m)
Info.MaxTerrainDist = 10000;            % Maximum Distance to look for shadows (m)
Info.NAngleDivisions = 36;              % Number of Divisions for the creation of the Solar Index
   
