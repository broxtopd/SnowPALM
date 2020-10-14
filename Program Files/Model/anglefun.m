function [alphamat] = anglefun(x,y,z,mask,Alphas,divmap,mindist,maxdst)
% Computes the horizon angles for each cell in each direction (used for the
% solar forcing and Sx maps)
%
% Inputs: x, y, and z are vectors of containing the x, y, and z coordinates of the
%           surface (elevation + vegetation height)
%         mask is a mask showing which locations to perform the analysis for
%         Alphas is a vector containing the angle increments in the analysis
%         divmap is a map showing how individual map tiles are broken up (so tiles
%           can be processed one at a time)
%         maxdst is the maximum search distance
% Outputs: alphamat - 2d vector of horizon angles om each direction
%
% Created by Patrick Broxton (broxtopd@email.arizona.edu)
% Updated June 2012

%#codegen
% % To improve performance, matlab coder can be used to create c++ code from
% % this function, and to compile as a mex function
% % To compile, type 
% To compile, type 
% cfg = coder.config;
% cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
% cfg.IntegrityChecks = false;
% cfg.ExtrinsicCalls = true;
% cfg.ResponsivenessChecks = true;
% codegen anglefun -config cfg -args {coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(double(0),[inf,inf]),coder.typeof(uint16(0),[inf,inf]),double(0),double(0)}
%
% % This line also needs to be uncommented if using matlab coder
% coder.extrinsic('disp', 'num2str');  

% Store result into an 8 bit matrix
alphamat = uint8(zeros([numel(x) numel(Alphas) - 1]));
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
    z2 = z(idx2);
    mask2 = mask(idx2);
         
    % ...as well as areas adjacent to the current map tile because they
    % affect the values of Sb indices in the current map tile
    x2a = x(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    y2a = y(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    z2a = z(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);

    % Now, for each node in these tiles, compute the horizon angles
    for id2 = 1:numel(idx2)
        % If this index is within the portion of the map where we want to
        % compute the angles
        if mask2(id2) == 1
            % Again, only search out to 'maxdst'
            x3 = (x2a - x2(id2));   y3 = (y2a - y2(id2)); z3 = (z2a - z2(id2));

            dist = (x3.^2 + y3.^2) .^ 0.5;
            locs = dist < mindist | dist > maxdst | z3 < 0;
            x3(locs) = [];      y3(locs) = [];    z3(locs) = [];
            tantheta = z3 ./ (x3.^2 + y3.^2) .^ 0.5;            

            x3(tantheta<0.01) = [];      y3(tantheta<0.01) = [];    z3(tantheta<0.01) = [];   tantheta(tantheta<0.01) = [];   
            tan_map = y3 ./ x3;
            tan_map_l = tan_map; tan_map_l(x3>=0) = NaN;
            tan_map_r = tan_map; tan_map_r(x3<0) = NaN;            
            
            for alphan = 1:numel(Alphas)-1
                % Find maximum and minimum (horizontal) angles for each direction
                max_alpha = Alphas(alphan+1);
                min_alpha = Alphas(alphan);
                cos_max_alpha = cos(max_alpha);
                cos_min_alpha = cos(min_alpha);
                tan_max_alpha = tan(max_alpha);
                tan_min_alpha = tan(min_alpha);
                pos_map = zeros(size(x3));    
            
                % Find the positions of all of the cells that are included in
                % that slice (depending on the sign of the cosine, decide the 
                % direction of the greater thans
                if cos_max_alpha >= 0 && cos_min_alpha >= 0
                    pos_map(tan_map_r >= tan_min_alpha & tan_map_r <= tan_max_alpha) = 1;
                elseif cos_max_alpha < 0 && cos_min_alpha >= 0
                    pos_map(tan_map_r >= tan_min_alpha | tan_map_l <= tan_max_alpha) = 1;
                elseif cos_max_alpha < 0 && cos_min_alpha < 0
                    pos_map(tan_map_l >= tan_min_alpha & tan_map_l <= tan_max_alpha) = 1;
                elseif cos_max_alpha >= 0 && cos_min_alpha < 0
                    pos_map(tan_map_l >= tan_min_alpha | tan_map_r <= tan_max_alpha) = 1;
                end
               
                % Find the maximum angle
                if ~isempty(tantheta(pos_map>0))
                    tmp = atan(max(tantheta(pos_map>0)));
                    alphamat(idx2(id2),alphan) = uint8(tmp*(255/(pi/2)));
                else
                    alphamat(idx2(id2),alphan) = uint8(0);
                end
            end
            c = c+1;
        else
            alphamat(idx2(id2),:) = uint8(0);
        end
    end
    % Display function progress
    if c/sum(mask) > progthresh
        disp(['Finding angles: ' num2str(round(c/sum(mask)*100)) '% complete']);
        progthresh = double(c/sum(mask)) + 0.05;
    end
end