function centers = filt_boundaries(img,centers, search_boxxR,search_boxyD,search_boxxL, search_boxyU)


    for k =1:size(centers,1)
        if ceil(centers(k,1)) > size(img,2)-search_boxxR-1  || ceil(centers(k,2)) > size(img,1)-search_boxyD-1 ...
            || search_boxxL + 1 > floor(centers(k,1)) || search_boxyU + 1 > floor(centers(k,2)) 
        centers(k,:)=nan;
        end
    end
    centers(any(isnan(centers), 2), :) = [];   

end

