function [minpeaks,maxpeaks, minmax_scaled] = horizprofile(xc,yc,img_normed,img_scaled,width)

    %for the number of peaks
    ROI_standardised  = sq_mask(medfilt2(img_normed, [5,5]),xc,yc,width);
    ROI_standardised =  ROI_standardised(2:end-1,2:end-1);
    LMin = islocalmin(mean(ROI_standardised,1), 'MinProminence',0.15);
    LMax = islocalmax(mean(ROI_standardised,1), 'MinProminence',0.15);
    minpeaks = sum(LMin);
    maxpeaks = sum(LMax);
    
    %for the profile of min max
    ROI_scaled  = sq_mask(imgaussfilt(img_scaled,10),xc,yc,width);
    ROI_scaled =  ROI_scaled(2:end-1,2:end-1);
    %square_scaled = mean(ROI_scaled,1);
    minmax_scaled = max(mean(ROI_scaled,1))-min(mean(ROI_scaled,1));


end
