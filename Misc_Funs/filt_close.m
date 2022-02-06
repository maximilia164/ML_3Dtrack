function centers = filt_close(centers,searchD)

Idx = rangesearch(centers,centers,searchD);

for m = 1:numel(Idx)
    if length(Idx{m}) >1
        centers(m,:) = nan;
    end
end

centers(any(isnan(centers), 2), :) = [];    

end
