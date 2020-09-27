function [canpos,wgtfun] = dstfun(x,y,vegcover,veght,mask,divmap,maxdst,dstpars)
% Finds each cells distance from the nearest canopy and builds a cone of
% influence out from the canopy (weighting larger trees more than smaller
% trees).
%
% Inputs: x and y - vectors of containing the x and y coordinates of the surface 
%           (elevation + vegetation height)
%         vegcover is a vectorized version of the vegetation cover map
%         veght is a vectorized version of the vegetation height map
%         divmap is a map showing how individual map tiles are broken up (so tiles
%           can be processed one at a time)
%         maxdst is the maximum search distance
%         dstpars is a vector containg a range of areas of influence of
%         radiation [m] (entire range is available for model execution, so
%         this can be adjusted using a model parameter)
% Outputs: canpos - vector showing the locations of the nearest canopy pixel
%          wgtfun - 2d vector with canopy cones generated for different
%          values of dstpars
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
% codegen dstfun -config cfg -args {coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(single(0),[inf,inf]),coder.typeof(uint16(0),[inf,inf]),double(0),coder.typeof(double(0),[inf,inf])}
%
% % This line also needs to be uncommented if using matlab coder
% coder.extrinsic('disp', 'num2str'); 

% Initialize the output variables
posmat = 1:numel(vegcover);
canpos = zeros(size(vegcover));
wgtfun = zeros([numel(dstpars) numel(vegcover)]);
% wgtfun = [];

idx = unique(divmap(:)');
progthresh=0;

c = 0;
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
    vegcover2 = vegcover(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    veght2 = veght(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    posmat2 = posmat(x>min(x(idx2))-maxdst & x<max(x(idx2))+maxdst & y>min(y(idx2))-maxdst & y<max(y(idx2))+maxdst);
    
    for id2 = 1:numel(idx2)
        % If this index is within the portion of the map where we want to
        % compute the angles
        if mask2(id2) == 1
            x3 = (x2a - x2(id2));   y3 = (y2a - y2(id2)); 
            dstmat = sqrt(x3.^2+y3.^2);
            dstmat_find = (max(dstmat(:)) - dstmat) .* vegcover2;
            Pos = posmat2(dstmat_find == max(dstmat_find(:)));
            % Find canopy positions
            canpos(idx2(id2)) = Pos(1);
            % Build cones out from the trees (e.g. a new surface where each
            % tree has an enlarged area of influence)
            for i = 1:numel(dstpars)
                wgtfun(i,idx2(id2)) = max(exp(-dstmat/dstpars(i)).*veght2);
            end
            c = c+1;
        else
            for i = 1:numel(dstpars)
                wgtfun(i,idx2(id2)) = uint8(0);
            end
        end
    end

    % Display function progress
    if c/sum(mask) > progthresh
        disp(['Finding distances to Vegetation: ' num2str(round(c/sum(mask)*100)) ' % complete']);
        progthresh = double(c/sum(mask)) + 0.05;
    end
end

