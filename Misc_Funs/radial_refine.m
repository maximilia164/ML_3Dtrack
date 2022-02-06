function out = radial_refine(im, pos, sz)

%Max Bailey 2021
%take the output of the IDL (Dufresne) and apply the radial centre finding
%method by Copyright 2011-2012, Raghuveer Parthasarathy, The University of Oregon
%PERFORMANCE: radial approach 100x faster than gaussian (not rate limiting)
%PURPOSE: Finding centres of BF particles more accurately, especially for
%out of focus (to be used with 3D diffraction tracking)
% NB: Does not handle dimers/trimers well

%im is the input image for tracking
%pos is the vectors of previously found by IDL 
%sz is the size of the mask around the coords to evaluate over
%%NOTE: RADIAL GRADIENT METHOD MADE MUCH BETTER BY RESCALING

im=double(im);
im = 255-im; %ALMOST ALWAYS BETTER TO USE A_INV
%inversion seems to help
if sz/2 == floor(sz/2)
warning('sz must be odd, like bpass');
end

if isempty(pos)
    warning('there were no positions inputted into cntrd. check your pkfnd theshold')
    out=[0 0 0 0 0 0];
    return;
end

%create mask - window around trial location over which to calculate the centroid

%round centroids for mask fitting
pos(:,1) = round(pos(:,1));
pos(:,2) = round(pos(:,2));

xmins = zeros(numel(pos(:,1)),1);
ymins = zeros(numel(pos(:,2)),1);

%initialise mask for each particle - on a frame basis on external for loop
for i = 1:numel(pos(:,1))
    
    im_temp = im; %keep loading image for each particle to crop
    
    %crop the image to fit the mask (computation time, important for radial)
    [sy, sx] = size(im_temp); %initial size - necessary each for loop?
    xmin = max(1, floor(pos(i,1)-sz)); %finds the relevant minimum pixel in x
    xmins(i) = xmin; %save the xmin
    xmax = min(sx, ceil(pos(i,1)+sz)); %finds the relevant maximum pixel in x
    ymin = max(1, floor(pos(i,2)-sz)); % ""y""
    ymins(i) = ymin;
    ymax = min(sy, ceil(pos(i,2)+sz)); % ""y""
    
    im_temp = im_temp(ymin:ymax, xmin:xmax); % trim boundaries of temporary image
    
    pos(i,1) = pos(i,1) - xmin + 1; %rescaling to new boundaries
    pos(i,2) = pos(i,2) - ymin + 1;
    
    %initialise mesh grid and mask
    [x, y] = meshgrid(1:size(im_temp,2), 1:size(im_temp,1)); 
    mask(:,:,i) = zeros(size(im_temp));
    
    %create mask with new boundaries and centres
    mask(:,:,i) = (x-pos(i,1)).^2 + (y-pos(i,2)).^2 < sz.^2;
    %multiply by temporary rescaled image
    mask(:,:,i) = mask(:,:,i).*double(im_temp);
    %DYNAMIC RESCALING AGAIN ON MASK (ENHANCE GRADIENTS FOR GRADIENT FINDING)
    mask(:,:,i) = 255.*(mask(:,:,i)-min(min(mask(:,:,i))))./(max(max(mask(:,:,i)))-min(min(mask(:,:,i))));
end

%remove particles/centroids within radius of mask from edge of image
[nr,nc]=size(im);
%remove all potential locations within distance sz from edges of image
%need to rescale with the xmins and ymins

ind=find((pos(:,2) + ymins -1) > (sz+2)/2 & (pos(:,2) + ymins -1) < nr-(sz+2)/2 & (pos(:,1) + xmins -1) ...
    > (sz+2)/2 & (pos(:,1) + xmins-1) < nc-(sz+2)/2);
%takes indices of particles which fulfill criteria of being within image
%including mask %y>size: for the 0 bound, y < edge-size - for top side
pos=pos(ind,:); %non normalised wrt mask offset
mask = mask(:,:,ind); 

% Np = length(pos); %number particles - doesn't work for 1 particle basis
Np = size(pos,1);

%loop through all of the candidate positions
pts=zeros(Np,4);

%condition that ensures particles in view
if Np == 0
    warning('No particles remaining after screening')
    out=[0 0 0 0 0 0];
    return;
end

if Np==1
    [pts(1), pts(2), pts(3), pts(4)] = radialcenter(mask(:,:));
else
    for j = 1:Np
      [pts(j,1), pts(j,2), pts(j,3), pts(j,4)] = radialcenter(mask(:,:,j));
%       [pts(j,1), pts(j,2), pts(j,3), pts(j,4)] = radialcenter(mask(:,:));
%     [pts(j,1), pts(j,2)] = radialcenter(mask(:,:,j));
    end
end

%renormalise to get back the shifted values for the mask
pts(:,1) = pts(:,1) + xmins -1;
pts(:,2) = pts(:,2) + ymins -1;

out = pts(:,1:2);
end

