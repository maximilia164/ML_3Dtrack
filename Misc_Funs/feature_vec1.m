function [M1]=feature_vec1(A)

%1st invariant moment of image

% First Moment
n20=cent_moment(2,0,A);
n02=cent_moment(0,2,A);
M1=n20+n02;
end
