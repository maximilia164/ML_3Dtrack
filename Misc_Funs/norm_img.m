function normed = norm_img(img)

    normed =  (img-mean(img(:)))/std(img(:));

end