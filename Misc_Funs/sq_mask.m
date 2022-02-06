function masked = sq_mask(img, xc, yc, r)
%mask of image - sq
%   M = meanDisk(IMG, XC, YC, R) returns the mean of IMG(Y,X) for all X and
%   Y such that the Euclidean distance of (X,Y) from (XC,YC) is less than
%   R. IMG must be 2-D, R must be positive, and some elements of IMG must
%   lie within the circle.
% This section is for efficiency only - avoids wasting computation time on
% pixels outside the bounding square
[sy,sx] = size(img);
xmin = max(1, floor(xc-r)); %role of sx/sy to take the maximum in the image
xmax = min(sx, ceil(xc+r));
ymin = max(1, floor(yc-r));
ymax = min(sy, ceil(yc+r));


img = img(ymin:ymax, xmin:xmax); % trim boundaries
xc = xc - xmin + 1;
yc = yc - ymin + 1;


% Make a square mask around image 
[x,y] = meshgrid(1:size(img,2), 1:size(img,1)); 
mask = (x-xc).^2 < r.^2 & (y-yc).^2 < r.^2; %mask with offset around centre
masked = mask.*double(img);
end
