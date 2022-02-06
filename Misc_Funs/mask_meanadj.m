function mask = mask_meanadj(img,set_mean)

    mean_val = mean(img(:));
    diff_mean = set_mean - mean_val;
    mask = img + diff_mean;

end
