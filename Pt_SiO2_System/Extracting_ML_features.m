clc; close all; clear all

%{
This script extracts the image features from the masks necessary for
training of "traditional" ML models. 
For Pt-SiO2 system, I adjusted the Z label
%}



%% Pipeline extracting image features from particle masks

umpp = 1/6.1538;
dp = 2.12 ;
diam = round(dp/umpp);
if mod(diam,2) == 0
    diam = diam + 1;
end

b_diam = 3*diam+10; %for the radial refine tracking step
mask_baseline = 145; %baseline mean of the particles masks set

% Parameters to evaluate image features over
tracking_diam = 2*diam+1;
horiz_vet = round(3*diam);
radial_sweep = 3*diam+1;

%{
Adjusting the Z-label
Load the radial profile extracted from an average of 5 particles, 
identified at the point where they become dark circles, labelled as 
approximately 5 um above the focal point 
%}
% load('*/5um_RadialProfile.mat')
Z_step = 0.1; %step size in um
Z = -20:Z_step:20; %Original Z-stack indices, before cropping etc.
Z_set = 5; 


%%
% pp=what('Particle_Masks');
listing=dir([pp.path '/*.mat']);

Dataset=[];
x_nc = {};
y_nc = {};
particle_numbers = [];

for p = 1:numel(listing) %iterate over the different Z stacks
    
%loading the particle positions + masks previously extracted
    particles = load(listing(p).name).particles; 
    particle_numbers(p) = length(particles);
    
    for i = 1:length(particles) %first loop through the particles in each Z-stack
        
        %First step is to re-index Z-values according to reference 5um label
        Z_labels = cell2mat(particles{i}(:,3)); %Extract the Z labels of the particle 
        Z_vals = Z(Z_labels);
        frame_num = size(particles{i},1);
        radial_check = zeros(frame_num,radial_sweep);
        x_nc{i} = zeros(frame_num,1);
        y_nc{i} = zeros(frame_num,1);

        for k = 1:frame_num
            a = particles{i}{k,5};
            xy = FindMiddle(a);
            refined2 = radial_refine(a,xy,tracking_diam); %improving the centre finding using custom tracking script
            x_nc{i}(k) = round(refined2(:,1));
            y_nc{i}(k) = round(refined2(:,2));
            a = rescale(a);
            [radial_check(k,:), ~] = radialprofile(round(refined2(:,1)),round(refined2(:,2)),a,radial_sweep); 
        end
        
        %Adjusting the Z-values
        Z_idx = Z_vals(min(least_squares_finding(radial_label,radial_check,10))); %the lowest frame closest to target
        diff_Z = Z_set - Z_idx;
        Z_vals = Z_vals + diff_Z; %adjusting as necessary

        for j = 1:frame_num %looping over Z-slices to extract features with corrected labels

            mask = particles{i}{j,5};
            mask = mask_meanadj(mask, mask_baseline); 
            scaled = rescale(mask);
            normed = norm_img(mask);

            xnc = x_nc{i}(j);
            ync = y_nc{i}(j);
            
            %Categorical Features and Horizontal Profile Information extraction
            [minpeaks,maxpeaks,minmax_scaled] = horizprofile(xnc,ync,normed,scaled,horiz_vet);
            
            %Radial Profile extraction
            mask = medfilt2(mask,[3,3])./255; 
            [radial, r_minmax] = radialprofile(xnc,ync,mask,radial_sweep);

            %Extracting Image Inertia
            mask_inertia =  sq_mask(mask,xnc,ync, diam); 
            mask_inertia = mask_inertia(2:end-1, 2:end-1);
            moments = feature_vec1(mask_inertia); 
            
            %All the features together
            data = [moments minmax_scaled minpeaks maxpeaks minpeaks+maxpeaks r_minmax radial Z_vals(j)];
            Dataset = [Dataset; data];

        end  

    end
end

