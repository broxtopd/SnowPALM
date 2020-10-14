function output_csv(Info,StartDate,EndDate,SubDomain,MVars)
% Read in .mat files containing tabular output data generated during the 
% model runs, combine them (they are generated for each tile separately),
% and write results to CSV files (both for hourly output, and for daily
% model outputs)
%
% Inputs: Info - Structure containing information about the model run
%         StartDate - Starting Date for the simulation
%         EndDate - Ending Date for the simulation
%         SubDomain - Model Domain name
%         MVars - Structure containing the names of the model variables on
%           the output tape
% No outputs, instead generates CSV files
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated April 2017

gen.HOUR = 3600;
gen.DAY = 86400;
gen.TS = gen.HOUR;                  % Model timestep = 1 hr
[Info,Cutout] = get_bounds(Info);   % Get the list of tiles
Info.NameIdentifier = SubDomain;
STStamp = datenum(StartDate);
ETStamp = datenum(EndDate);
TSVec = STStamp:gen.TS/gen.DAY:ETStamp+1-(gen.TS/gen.DAY);
TSVec_Day = STStamp:ETStamp;

% Make the output directory if it does not exist
if ~exist([Info.DisplayDir filesep 'Tabular' filesep Info.NameIdentifier],'file')
    mkdir([Info.DisplayDir filesep 'Tabular' filesep Info.NameIdentifier]);
end

% Load output from the first available tile to get the list of locations
% that have data
exec = 0;
Info2 = Info;
for i = 1:numel(Cutout)
    Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(i)];
    if exist([Info.DisplayDir filesep 'Tabular' filesep Info2.NameIdentifier '.mat'],'file')
        load([Info.DisplayDir filesep 'Tabular' filesep Info2.NameIdentifier '.mat']);
        exec = 1;
        break
    end
end

if exec == 1
    % Loop through the locations
    for i = 1:numel(TabularData)
        OutVars = [];
        % Record the list of variable names on the output tape
        for v = 1:numel(MVars)
            OutVars{v} = MVars(v).Name;
        end
        OutMat = [];
        TotalArea = 0;
        % Loop through the tiles...
        for h=1:numel(Cutout)
            Info2 = Info;
            Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(h)];
            Info2.Cutout = Cutout(h);
            % Load the .mat file containing the tabular data ..., if a POI
            % (really area of interest) falls on several tiles, keep track
            % of the weighted average of model output variables on all tiles
            if exist([Info.DisplayDir filesep 'Tabular' filesep Info2.NameIdentifier '.mat'],'file')
                load([Info.DisplayDir filesep 'Tabular' filesep Info2.NameIdentifier '.mat']);
                % If a POI is on a particular tile (i.e. the POI's area on the
                % tile is > 0), data from the tile contributes to the weighted
                % average
                if sum(TabularData(i).area) > 0
                    TotalArea = TotalArea + TabularData(i).area;
                    if isempty(OutMat)
                        for v = 1:numel(MVars)
                            eval(['OutMat(:,v) = TabularData(i).' MVars(v).Name ' * TabularData(i).area;']);  
                        end
                    else
                        for v = 1:numel(MVars)
                            eval(['OutMat(:,v) = OutMat(:,v) + TabularData(i).' MVars(v).Name ' * TabularData(i).area;']);  
                        end
                    end
                end
            end
        end
        % At this point, we have a weighted sum

        if ~isempty(OutMat)
            % Convert to weighted average
            OutMat = OutMat ./ TotalArea;
            % Compute daily values from hourly values
            for u = 1:numel(OutMat(:,1))/24
                OutMat_Day(u,:) = mean(OutMat((u-1)*24+1:u*24,:),1);
            end

            % Write the output CSV files
            fid = fopen([Info.DisplayDir filesep 'Tabular' filesep Info.NameIdentifier filesep char(TabularData(i).Name) '_ts.csv'],'w');
            fprintf(fid,'Timestep,');
            for t = 1:numel(OutVars)
                fprintf(fid,OutVars{t});
                fprintf(fid,',');
            end
            fprintf(fid,'\n');
            for u = 1:numel(OutMat(:,1))
                fprintf(fid,'%s,',datestr(TSVec(u)));
                fprintf(fid,'%0.5f,',OutMat(u,:));
                fprintf(fid,'\n');
            end
            fclose(fid);

            fid = fopen([Info.DisplayDir filesep 'Tabular' filesep Info.NameIdentifier filesep char(TabularData(i).Name) '_daily.csv'],'w');
            fprintf(fid,'Date,');
            for t = 1:numel(OutVars)
                fprintf(fid,OutVars{t});
                fprintf(fid,',');
            end
            fprintf(fid,'\n');
            for u = 1:numel(OutMat_Day(:,1))
                fprintf(fid,'%s,',datestr(TSVec_Day(u)));
                fprintf(fid,'%0.5f,',OutMat_Day(u,:));
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
    end

    % Delete Tabular Data .mat files that were generated during the model runs
    for h=1:numel(Cutout)
        Info2.NameIdentifier = [Info.NameIdentifier '_' num2str(h)];
        if exist([Info.DisplayDir filesep 'Tabular' filesep Info2.NameIdentifier '.mat'],'file')
            delete([Info.DisplayDir filesep 'Tabular' filesep Info2.NameIdentifier '.mat']);
        end
    end
end