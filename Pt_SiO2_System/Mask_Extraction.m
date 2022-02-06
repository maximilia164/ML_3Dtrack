clc; close all; clear all


%{
This is the first script of tracking. It extracts out the particle masks
from the Z-stacks for training the ML models (Pt/SiO2 system). It works
because across Hough Circle transform, I used the same Z-stack range. For
the TiO2-SiO2 system, I adjust the range, and therefore a slightly modified
script was used

%}
%% Loading image stacks
% pp=what('MBN19_T1F0_Zstack');
listing=dir([pp.path '/*.tif']);
img ={}; 
img_orig = {};

for j=1:numel(listing) %for loop through all the sims
    img{j} = double(imread(listing(j).name));
    img_orig{j} = img{j};
    img{j} = 255.*(img{j}-min(min(img{j})))./(max(max(img{j}))-min(min(img{j})));
end
%% Image Tracking Testing

%parameters of the microscope and particles
umpp = 1/6.1538;
dp = 2.12 ;
diam = round(dp/umpp);
if mod(diam,2) == 0
    diam = diam + 1;
end

%parameters for Hough Circle
sensitivity=0.9;
polarity='bright';
method1='twostage';

%filtering par
search_box = 3*diam; 
searchd = ceil(sqrt(2*search_box^2)); %size to sweep based on masks used - make sure only one particle in each mask

%check the tracking
a = img{250};
a_gauss = medfilt2(a,[5,5]);
a_inv = 255-a_gauss;
% a_inv = rescale(a_inv);
b_hough = bpass(a_inv,1,diam*3);
[centersBright, radiiBright] = imfindcircles(b_hough,[floor(diam/2) diam*2],'ObjectPolarity',polarity,'Sensitivity',sensitivity,'Method',method1 );
centers = filt_close(centersBright,searchd);
figure()
colormap('gray'), imagesc(a);
hold on
plot(centers(:,1),centers(:,2),'r*', 'DisplayName', 'original')
axis square;

%% Step 1: Particle centre finding
pxpy=[];
masks = {};

start_vid = 151; %where particle tracking starts to be consistent (across Z)
end_vid   = 376; %where particle tracking stops to be consistent

umpp = 1/6.1538;
dp = 2.12 ;
diam = round(dp/umpp);
if mod(diam,2) == 0
    diam = diam + 1;
end
sensitivity=0.9;
polarity='bright';
method1='twostage';

pxpy=[];
search_box = 3*diam; %based on horizvet - also good for hough filt bpass
b_diam = 2*diam; %seaching diameter of the hough transform etc.
searchd = ceil(sqrt(2*search_box^2))+10; %size to sweep based on masks used - make it larger in case centres are not exactly aligned

for i = start_vid:end_vid %load across the frames found to have valid centres found
    
    a = img{i}; %load image from the stack
    a_filt = medfilt2(a,[5,5]); %median filt works better for edge preservation 
    a_inv = 255-a_filt; 
    b_hough = bpass(a_inv,1,b_diam); %bandpass applied for better centre finding - using filtered images

    %Hough Centre Finding
    [centers, ~] = imfindcircles(b_hough,[floor(diam/2) b_diam],'ObjectPolarity',polarity,'Sensitivity',sensitivity,'Method',method1 );
    centers(:,3) = i; %assigning the frame

    % filtering out the particles too close to boundaries or to each other
    % important because full scattering rings required without interference
    centers = filt_close(centers,searchd+10); %filters out the centres too close to each other  
    centers = filt_boundaries(a,centers,search_box,search_box,search_box,search_box); %filters out particles too close to edge for box needed
    
    pxpy = [pxpy; centers];

end

%% Linking trajectories and separating particles

param.mem = 5; %since Z stack - unlikely particles randomly switch
param.dim = 2;
param.good = 150; %delete particles contributing to less than 150 frames (removes noisy)
param.quiet=1;
fps = 10;
max_d = ceil(5/fps/umpp); %5ums-1, divided by fps, divided by umpp for pixel value

tracked = track(pxpy,max_d, param);
pids = unique(tracked(:,4));
particles = {};

for i = 1:numel(pids)
    particles{i} = tracked(any(tracked(:,4)==pids(i),2),:); %collect particles into cells
end

%% checking image tracking - make sure Z stacks ok
figure()
a = img{201};
colormap('gray'), imagesc(a);
hold on
for i = 1:numel(particles)
    plot(particles{i}(:,1), particles{i}(:,2),'linewidth',1.5, 'LineStyle','-')
end
axis square;

%% adding mask for each particle as 5th column

%make each particle it's own cell
for i = 1:length(particles)
    particles{i} = num2cell(particles{i}); 
end

%extracting masks for each particle
for i = 1:length(particles) %loop first over particles
    
    for j = 1:size(particles{i},1) %looping over frames
        
        %getting details out of cell
        xc = round(particles{i}{j,1}); %xcoord
        yc = round(particles{i}{j,2}); %ycoord
        frame = particles{i}{j,3}; %frame number
        
        %mask extraction - tied to each particle
        a_orig = img_orig{frame};
        mask = sq_mask(a_orig,xc,yc,searchd); %getting mask out from relevant frame of particle from original image
        mask =  mask(2:end-1,2:end-1); %cropping the mask - removing zeros
        particles{i}{j,5} = mask;
       
    end    
        
end

%% sanity check on number of particles over frames
histogram(pxpy(:,3), length(img))
% a = img{354};
% figure()
% colormap('gray'), imagesc(a);
