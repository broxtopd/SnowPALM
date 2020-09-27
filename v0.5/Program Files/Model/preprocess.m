function preprocess(Info,SubDomain)
% Function to, if necessary, precompute the spatial indexes, and assemble
% the Sx map.  This only needs to occur if the Sx map has not been
% assembled yet (first, it will be created in tiles), and then it has to be
% mosaiced before Sb calculations can be done.
%
% Inputs: Info - Structure containing information about the model run
%         SubDomain - Model Domain name
% No outputs, instead generates a Sx map (as a GIS file)
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated July 2017

[Info] = program_pars(SubDomain,Info);   % Get spatial parameters
PWD = pwd;
cd(Info.ProgramFilesDir)
set_paths(Info.ProgramFilesDir,'add')
cd(PWD)

[modelpars] = model_pars;
Info.NameIdentifier = SubDomain;
    
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

% Name of the Sx File
sxfile = [Info.SpatialDataDir filesep SubDomain '_Sx' filesep 'Sx_' num2str(mindir) '_' num2str(maxdir) '_' num2str(modelpars.Sx_Exponent) '_' num2str(modelpars.Sx_Rescale) '.tif'];
% For this step, make sure that we are not computing Sb
ComputeSbIndex = Info.ComputeSbIndex;
Info.ComputeSbIndex = 0;

% Only create the Sx file if it is to be computed for this model subdomain
% and it does not exist
if Info.ComputeSxIndex && strcmp(SubDomain,Info.SxSource) && ~exist(sxfile,'file')
    [Info,Cutout] = get_bounds(Info);
    % Run the process_spatial program (this will make sure that all
    % indexes, including Sx, but not Sb, are created)
    % Built in Parallelization
    if Info.RunParallel == 1 && numel(Cutout) > 1 && Info.MaxNumWorkers > 1          
        try 
            parpool(min(Info.MaxNumWorkers,numel(Cutout))); 
        end
        parfor i=1:numel(Cutout)
            Info2 = Info;
            Info2.NameIdentifier_orig = Info.NameIdentifier;
            Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
            Info2.Cutout = Cutout(i);
            [SpatialData] = process_spatial(Info2,modelpars);
            IdentifierList(i).Identifier = Info2.NameIdentifier;
        end
        pause(2)
    % Serial Mode 
    else        
        for i=1:numel(Cutout)
            Info2 = Info;
            Info2.NameIdentifier_orig = Info.NameIdentifier;
            Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
            Info2.Cutout = Cutout(i);
            [SpatialData] = process_spatial(Info2,modelpars);
            IdentifierList(i).Identifier = Info2.NameIdentifier;
        end
        pause(2)
    end

    combine_sx_maps(Info,IdentifierList,modelpars);
end

Info.ComputeSbIndex = ComputeSbIndex;
