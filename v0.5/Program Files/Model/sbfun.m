function [smat] = sbfun(x,y,mask,alphamat,Alphas,divmap,maxdst)
% Computes the Sb index that is used in the Winstrel wind model
%
% Inputs: x and y - vectors of containing the x and y coordinates of the surface 
%           (elevation + vegetation height)
%         mask is a mask showing which locations to perform the analysis for
%         alphamat is the Sx "surface" that is used in the calculation of Sb values
%         Alphas is a vector containing the angle increments in the analysis
%         divmap is a map showing how individual map tiles are broken up (so tiles
%           can be processed one at a time)
%         maxdst is the maximum search distance
% Outputs: smat - vector of Sx values for a given wind direction
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated 6/2012

%#codegen
% % To improve performance, matlab coder can be used to create c++ code from
% % this function, and to compile as a mex function
% % To compile, type 
% cfg = coder.config;
% cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
% cfg.IntegrityChecks = false;
% cfg.ExtrinsicCalls = true;
% cfg.ResponsivenessChecks = true;
% codegen sbfun -config cfg -args {coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(uint8(0),[inf,inf]),coder.typeof(double(0),[inf,inf]),coder.typeof(uint16(0),[inf,inf]),double(0)}
%
% % This line also needs to be uncommented if using matlab coder
% coder.extrinsic('disp', 'num2str'); 

% Store result into an 8 bit matrix
smat = int8(zeros([numel(x) numel(Alphas) - 1]));
% Get unique indexes for each of the map tiles
idx = unique(divmap(:)');
Alphas = Alphas - pi;

progthresh=0;       % For displaying the function progress

c = 0;
% Loop through all of the map tiles ...
for id = 1:numel(idx)
    % ...select the appropriate portions of all of the necessry variables
    idx2 = find(idx(id) == divmap);
    x2 = x(idx2);
    y2 = y(idx2);
    mask2 = mask(idx2);
    % ...as well as areas adjacent to the current map tile because they
    % affect the values of Sb indices in the current map tile
    x2a = x(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    y2a = y(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    alphamat2a = alphamat(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst,:);
    
    % Now, for each node in these tiles, compute the Sx index
    for id2 = 1:numel(idx2)
        % If this index is within the portion of the map where we want to
        % compute Sb indices
        if mask2(id2) == 1
            % Again, only search out to 'maxdst'
            x3 = (x2a - x2(id2));   y3 = (y2a - y2(id2));
            locs = abs(x3) > maxdst | abs(y3) > maxdst;
            x3(locs) = [];      y3(locs) = [];

            tan_map = y3 ./ x3;
            tan_map_l = tan_map; tan_map_l(x3>=0) = NaN;
            tan_map_r = tan_map; tan_map_r(x3<0) = NaN;
            for alphan = 1:numel(Alphas)-1
                Sx = single(alphamat2a(~locs,alphan))/(255/(pi/2));
                SxVal = single(alphamat(idx2(id2),alphan))/(255/(pi/2));
                
                % Find maximum and minimum (horizontal) angles for each direction
                max_alpha = Alphas(alphan+1);
                min_alpha = Alphas(alphan);
                cos_max_alpha = cos(max_alpha);
                cos_min_alpha = cos(min_alpha);
                tan_max_alpha = tan(max_alpha);
                tan_min_alpha = tan(min_alpha);
                
                % Find the positions of all of the cells that are included in
                % that slice (depending on the sign of the cosine, decide the 
                % direction of the greater thans
                pos_map = zeros(size(x3));
                if cos_max_alpha >= 0 && cos_min_alpha >= 0
                    pos_map(tan_map_r >= tan_min_alpha & tan_map_r <= tan_max_alpha) = 1;
                elseif cos_max_alpha < 0 && cos_min_alpha >= 0
                    pos_map(tan_map_r >= tan_min_alpha | tan_map_l <= tan_max_alpha) = 1;
                elseif cos_max_alpha < 0 && cos_min_alpha < 0
                    pos_map(tan_map_l >= tan_min_alpha & tan_map_l <= tan_max_alpha) = 1;
                elseif cos_max_alpha >= 0 && cos_min_alpha < 0
                    pos_map(tan_map_l >= tan_min_alpha | tan_map_r <= tan_max_alpha) = 1;
                end

                % Find the difference between the current cell Sx and those
                % in the search domain
                tmp = SxVal - Sx(pos_map > 0);
                
                if ~isempty(tmp)
                    smat(idx2(id2),alphan) = int8(sum(tmp(:))/numel(tmp)*(127/(pi/2)));
                else
                    smat(idx2(id2),alphan) = int8(0);
                end
            end
            c = c+1;
        else
            smat(idx2(id2),:) = int8(0);
        end
    end

    % Display function progress
    if c/sum(mask) > progthresh
        disp(['Finding Sb index: ' num2str(round(c/sum(mask)*100)) ' % complete']);
        progthresh = double(c/sum(mask)) + 0.05;
    end
end