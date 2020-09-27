function [ulx,uly,lrx,lry,outdir,demfile] = domains(Domain)
% Define forcing data domains

%% Begin User Input

outdir = pwd;                                   % Local Directory where forcing data should be stored
demfile = 'US_DEM.tif';                      	% DEM to use along with the forcing data (will be reprojected)

% Forcing domains (Download areas that are larger than, or encompass
% multiple individual SnowPALM domains

Domains(1).name = 'Valles';
Domains(1).ulx = -107;
Domains(1).uly = 36.5;
Domains(1).lrx = -106;
Domains(1).lry = 35.5;

%% End User Input

num = 0;
for i = 1:numel(Domains)
    if strcmp(Domains(i).name,Domain)
        ulx = num2str(Domains(i).ulx);
        lrx = num2str(Domains(i).lrx);
        uly = num2str(Domains(i).uly);
        lry = num2str(Domains(i).lry);
        num = i;
    end
end

if num == 0
    error(['No Forcing Domain called ' Domain]);
end
    
