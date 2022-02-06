function idx = least_squares_finding(vec, mat, lowest)

[~,idx]=mink(sum((mat-vec).^2,2), lowest);

end
