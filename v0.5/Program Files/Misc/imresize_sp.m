function B = imresize_sp(A,scale,method)
% Closely matches the imresize function from the image processing toolbox
% Created by Patrick Broxton

if nargin < 3
    method = 'linear';
end

if strcmp(method,'bilinear')
    method = 'linear';
end

% Define X and Y matricies for the original data
sz_old = size(A);
X_old = repmat(1:sz_old(2),[sz_old(1) 1]);
Y_old = flipud(repmat((1:sz_old(1))',[1 sz_old(2)]));
X_old_nd = flipud(X_old)';
Y_old_nd = flipud(Y_old)';
A_nd = double(flipud(A)');

% Define X and Y matricies for the new data
if numel(scale) == 1
    sz_new = size(A) * scale;
    scale_x = scale;
    scale_y = scale;
elseif numel(scale) == 2
    sz_new = scale;
    scale = scale ./ size(A);
    scale_x = scale(2);
    scale_y = scale(1);
end
X_new = repmat(linspace(1-1/round(scale_x),sz_old(2)+1/round(scale_x),sz_new(2)),[sz_new(1) 1]);
Y_new = flipud(repmat((linspace(1-1/round(scale_y),sz_old(1)+1/round(scale_y),sz_new(1)))',[1 sz_new(2)]));
X_new_nd = flipud(X_new)';
Y_new_nd = flipud(Y_new)';

% Do the interpolation
interp = griddedInterpolant(X_old_nd,Y_old_nd,A_nd,method);
B_nd = interp(X_new_nd,Y_new_nd);
B = flipud(B_nd');
B = cast(B,'like',A);