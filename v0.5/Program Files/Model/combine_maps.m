function combine_maps(gen,Info,IdentifierList,FVars,MVars)
% Combines all of the map tiles from the model tile data created by the 
% display_output program
%
% Inputs: gen - general constants used by the model
%         Info - Structure containing information about the model run
%         IdentifierList - List of map tiles to combine
%         FVars - Structure containing the names of forcing variables on
%           the output tape
%         MVars - Structure containing the names of the model variables on
%           the output tape
% No outputs, instead generates GIS files
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)

STStamp = datenum(Info.StartDate);
ETStamp = datenum(Info.EndDate);

% Combines mosaics of General Terrain Maps
if Info.DEM
    combine_vrt_tif(Info,IdentifierList,'elev',[]); 
end
if Info.Slope
    combine_vrt_tif(Info,IdentifierList,'slope',[])
end
if Info.Northness
    combine_vrt_tif(Info,IdentifierList,'northness',[])
end
if Info.VegHT
    combine_vrt_tif(Info,IdentifierList,'VegHT',[])
end
if Info.Cover
    combine_vrt_tif(Info,IdentifierList,'cdensity',[]);
end
if Info.SkyView
    combine_vrt_tif(Info,IdentifierList,'skyview',[]);
end

% Solar forcing maps
if Info.SFI_BareEarth
    if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep 'sfi'],'file')
        mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep 'sfi'])
    end
    fname = ['sfi' filesep datestr(STStamp,'yyyymmdd') '_' datestr(ETStamp,'yyyymmdd')];
    combine_vrt_tif(Info,IdentifierList,fname,[]);
end

if Info.SFI_Veg
    if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep 'sfi_veg'],'file')
        mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep 'sfi_veg'])
    end
    fname = ['sfi_veg' filesep datestr(STStamp,'yyyymmdd') '_' datestr(ETStamp,'yyyymmdd')];
    combine_vrt_tif(Info,IdentifierList,fname,[]);
end

% Wind indices
if Info.Sx
    fname = 'Sx';
    combine_vrt_tif(Info,IdentifierList,fname,[]);
end

if Info.Sb
    fname = 'Sb';
    combine_vrt_tif(Info,IdentifierList,fname,[]);
end

% Loop through the days and make output maps for the forcing data and model
% outputs

parfor TS = STStamp:ETStamp     % V0.6: Change this to use parfor 
    for j = 1:gen.DAY/gen.TS
        TS_sub = TS+((j-1)/(gen.DAY/gen.TS));
        for v = 1:numel(FVars)
            if FVars(v).map && strcmp(Info.outputType,'ts')
                if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep FVars(v).Name '_ts'],'file')
                    mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep FVars(v).Name '_ts'])
                end
                combine_vrt_tif(Info,IdentifierList,MVars(v).Name,TS_sub); 
            end
        end
        for v = 1:numel(MVars)
            if MVars(v).map && strcmp(Info.outputType,'ts')
                if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep MVars(v).Name '_ts'],'file')
                    mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep MVars(v).Name '_ts'])
                end
                combine_vrt_tif(Info,IdentifierList,MVars(v).Name,TS_sub); 
            end
        end
    end

    for v = 1:numel(FVars)
        if FVars(v).map && strcmp(Info.outputType,'day')
            if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep FVars(v).Name '_day'],'file')
                mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep FVars(v).Name '_day'])
            end
            combine_vrt_tif(Info,IdentifierList,FVars(v).Name,TS_sub); 
        end
    end
    for v = 1:numel(MVars)
        if MVars(v).map && strcmp(Info.outputType,'day')
            if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep MVars(v).Name '_day'],'file')
                mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep MVars(v).Name '_day'])
            end
            combine_vrt_tif(Info,IdentifierList,MVars(v).Name,TS_sub); 
        end
    end       
end

% Finally, delete the directory containing the individual map tiles
for i = 1:numel(IdentifierList)
    if exist([Info.DisplayDir filesep 'GIS' filesep IdentifierList(i).Identifier],'file')
        rmdir([Info.DisplayDir filesep 'GIS' filesep IdentifierList(i).Identifier],'s');
    end
end


function combine_vrt_tif(Info,IdentifierList,Name,TS_sub)
% Builds vrt's out of all of the tiles generated for each model tile, and
% then converts to output GIS format (currently geotiff, but could be
% modified to be any / all gdal supported formats)

% Create output GIS folder if it does not exist
if ~exist([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier],'file')
    mkdir([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier]);
end

if ~isempty(TS_sub)
    eval(['typename = [''' Name '_day'' filesep datestr(TS_sub,''yyyymmdd'')];']);
else
    typename = Name;
end

if strcmp(Info.outputFType,'tif')
    disp(['Combining tiles for ' Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.tif']);
    % Create a list of all of the files to put together
    flist = [tempname '.txt'];
    fid = fopen(flist,'w');
    for i = 1:numel(IdentifierList)
      	addr = [Info.DisplayDir filesep 'GIS' filesep IdentifierList(i).Identifier filesep typename '.tif'];
        if exist(addr,'file')
            fprintf(fid,[strrep(addr,'\','/') '\n']);
        end
    end
    fclose(fid);

    % Create the tiled mosaic and convert to geotiff
    evalc(['!gdalbuildvrt -hidenodata -vrtnodata nan -input_file_list "' flist '" "' Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.vrt"']);
    evalc(['!gdal_translate -co COMPRESS=DEFLATE -co NUM_THREADS=ALL_CPUS "' Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.vrt" "' Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.tif"']);
    % The following code adds overview tiles (for fast GIS display of large datasets)
    if Info.CreateOverviewTiles
        evalc(['!gdaladdo "' Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.tif" 4 16 64']);
    end
    % Delete all of the temporary files
    delete([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.vrt']);
    delete(flist);
elseif strcmp(Info.outputFType,'mat')
    disp(['Combining tiles for ' Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.mat']);
    ny = round((Info.NSWE(1)-Info.NSWE(2))./Info.NewResolution);
    nx = round((Info.NSWE(4)-Info.NSWE(3))./Info.NewResolution);
    NSWE = Info.NSWE;
    Resolution = Info.NewResolution;
    Data = NaN(ny,nx);
    for i = 1:numel(IdentifierList)
        min_x = round(Info.Cutout(i).min_x0-Info.NSWE(3)./Info.NewResolution)+1;
        max_x = round(Info.Cutout(i).max_x0-Info.NSWE(3)./Info.NewResolution);
        min_y = round(Info.Cutout(i).min_y0-Info.NSWE(2)./Info.NewResolution)+1;
        max_y = round(Info.Cutout(i).max_y0-Info.NSWE(2)./Info.NewResolution);
        load([Info.DisplayDir filesep 'GIS' filesep IdentifierList(i).Identifier filesep typename '.mat']);
        Data(min_y:max_y,min_x:max_x) = flipud(val);
    end
    Data = flipud(Data);
    save([Info.DisplayDir filesep 'GIS' filesep Info.NameIdentifier filesep typename '.mat'],'-v7.3','Data','NSWE','Resolution');
end
   