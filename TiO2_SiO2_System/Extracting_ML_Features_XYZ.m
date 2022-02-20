clc; close all; clear all
%{
This script extracts out the ML features necessary to input into the
trained models from the particle masks previously found
%}
%%
umpp = 1/6.1538;
dp = 2.12 ;
diam = round(dp/umpp);
if mod(diam,2) == 0
    diam = diam + 1;
end
b_diam = 3*diam+10;
mask_baseline = 145;

tracking_diam = 2*diam+1;
horiz_vet = round(3*diam);
radial_sweep = 3*diam+1;



%%
particles_new = {};
for i = 1:length(particles) %first loop through the particles
    
    particles_new{i}(:,1) = cell2mat(particles{i}(:,1)); %setting up x in new particle cell
    particles_new{i}(:,2) = cell2mat(particles{i}(:,2)); %setting up y in new particle cell
    
    
    for j = 1:size(particles{i},1) %looping over frames
        
        
        %getting details out of cell
        xc = round(particles{i}{j,1}); %xcoord
        yc = round(particles{i}{j,2}); %ycoord
        mask = particles{i}{j,5};
        
        %getting out the masks 
        mask = mask_meanadj(mask, mask_baseline); %adjust the image mean to be setpoint

        %retrack centres better from the mask itself 
        xy = FindMiddle(mask);
        refined2 = radial_refine(mask,xy,tracking_diam);

        xnc = round(refined2(:,1));
        ync = round(refined2(:,2));
        
        scaled = rescale(mask);
        normed = norm_img(mask);

        %horizontal profile %remove square_scaled
        [minpeaks,maxpeaks,minmax_scaled] = horizprofile(xnc,ync,normed,scaled,horiz_vet);

        mask = medfilt2(mask,[3,3])./255; %rescaling the image after median filter blur
        [radial, r_minmax] = radialprofile(xnc,ync,mask,radial_sweep);

        %extracting 1st moment 
        mask_inertia =  sq_mask(mask,xnc,ync, diam); 
        mask_inertia = mask_inertia(2:end-1, 2:end-1);
        moments = feature_vec1(mask_inertia);
        data = [moments minmax_scaled minpeaks maxpeaks minpeaks+maxpeaks r_minmax radial];
   
        particles_new{i}(j,3:3+length(data)-1) = data;
 
       
    end  
    
end
