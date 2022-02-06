function [radial, r_minmax] = radialprofile(xc,yc,img,r_vet)

    radial = zeros(1,r_vet);
    
    for i = 1:r_vet
        radial(i) = meanDisk(img, xc, yc, i);
    end
    
    min_r = min(radial);
    max_r = max(radial);
    r_minmax = max_r-min_r;
    
end
