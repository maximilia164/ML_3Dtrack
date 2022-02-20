clc; close all; clear all
%%

%{
This script extracts out the masks of swimming particles, for processing as
necessary in the next steps. For XYZ tracking, stick to decision tree ML
approach for now (since models trained already)

%}


% pp=what('MBN19_T1F0_H2O2_3v_1');
listing=dir([pp.path '/*.tif']);

% first - load the images and find the means/differences
img ={}; 
img_orig = {};

%KEEP IT SIMPLE STUPID APPROACH
for j=1:numel(listing) %for loop through all the sims
    img{j} = double(imread(listing(j).name));
    img_orig{j} = img{j};
    img{j} = 255.*(img{j}-min(min(img{j})))./(max(max(img{j}))-min(min(img{j})));
end

%% Step 1: Particle centre finding
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
search_box = 3*diam+1; %based on horizvet - also good for hough filt bpass
searchd = ceil(sqrt(2*search_box^2)); 

for i = 1:numel(img)%_new) %img_new
    a = img{i}; %load image from set
    a_orig = img_orig{i};
    a_filt = medfilt2(a,[5,5]); %median filt works better for edge detection
    a_inv = 255-a_filt; 
    b_hough = bpass(a_inv,1,search_box);

    [centers, ~] = imfindcircles(b_hough,[floor(diam/2) diam*2],'ObjectPolarity',polarity,'Sensitivity',sensitivity,'Method',method1 );
    centers(:,3) = i;

    centers = filt_close(centers,searchd+10); %filters out the centres too close to each other  
    %!!NEED TO FILTER OUT THE BOUNDARIES LARGER - MESSES UP TRACKING 
    centers = filt_boundaries(a,centers,search_box+10,search_box+10,search_box+10,search_box+10); %filters out particles too close to edge for box needed
    pxpy = [pxpy; centers];

end


%% Linking trajectories and separating particles

param.mem = 2;
param.dim = 2;
param.good = 100; %delete with less than 50 frames (5s data)
param.quiet=1;
fps = 10;
max_d = ceil(40/fps/umpp); %20ums-1, divided by fps, divided by umpp for pixel value
% param = [2,2,5,1];
tracked = track(pxpy,max_d, param);
pids = unique(tracked(:,4));
particles = {};
for i = 1:numel(pids)
    
    particles{i} = tracked(any(tracked(:,4)==pids(i),2),:); %collect particles into cells
    %particles{i} = num2cell(particles{i});   
end

%% checking image tracking
figure()
a = img{1};
colormap('gray'), imagesc(a);
hold on
for i = 13%1:numel(particles)
    caption = num2str(i);
    plot(particles{i}(:,1), particles{i}(:,2),'linewidth',1.5, 'LineStyle','--', 'DisplayName',caption)
end
legend()

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
        %mask extraction
        a_orig = img_orig{frame};
        mask = sq_mask(a_orig,xc,yc,searchd); %getting mask out from relevant frame of particle from original image
        mask =  mask(2:end-1,2:end-1); %cropping the mask
        particles{i}{j,5} = mask;
       
    end    
        
end

%% sanity check - distribution of particle positions
histogram(pxpy(:,3), length(img))
%it is ok that there is a drop now - because at a certain point the
%particles are lost (high Z) - if they cannot be tracked properly, then
%that is a consequence of the 3D motion - better than fake values
%% sanity check of the mask
figure()
colormap('gray'), imagesc(particles{1}{1,5});
% hold on
% plot(centersBright(:,1),centersBright(:,2),'b*', 'DisplayName', 'original')

%% checking image tracking after changing formatting USE THIS TO CHECK TRAJS/Z
figure()
a = img{1};
colormap('gray'), imagesc(a);
hold on
for i = 12%48:48%numel(particles_new)
    caption = num2str(i);
    plot(particles_new{i}(:,1), particles_new{i}(:,2),'linewidth',1.5, 'LineStyle','--', 'DisplayName',caption)
end
legend()



