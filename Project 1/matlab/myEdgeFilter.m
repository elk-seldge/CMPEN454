function [img1] = myEdgeFilter(img0, sigma)
%Your implemention
h_size = 2 * ceil(3 * sigma) + 1;
gaussian_filter = fspecial("gaussian", h_size, sigma);
x_sobel = [-1 0 1; 
           -2 0 2; 
           -1 0 1];
y_sobel = [-1 -2 -1; 
            0 0 0; 
            1 2 1];
% apply the 'myImageFilter'
img0Smthd = myImageFilter(img0, gaussian_filter);
% cnvlv smthd img0 with sblX and sblY, get imgX and imgY df/dx df/dy
imgx = myImageFilter(img0Smthd, x_sobel);
imgy = myImageFilter(img0Smthd, y_sobel);
% convert rad to degree for convinient view
angl = atan2(imgy, imgx).* (180/pi);
% get magnitude from imgX and imgY
magnitude = sqrt(imgx.^2 + imgy.^2);
im_rec = img0; % nms perform on the replication of img0 then pass to img1

% NMS get rid of non-local-maximum in gradient neighbours
% angles: 0 45 90 135
angl(angl<0) = angl(angl<0) + 180; % set all angles to [0, 180] 
[row, col] = size(angl); % angle() calculate the angle and mag
for i = 2: row - 1
    for j = 2: col - 1
        if (angl(i, j) < 22.5) || (157.5 < angl(i, j)) % 0 degree ([0,45/2] union [(180-135)/2, 180])
            angl(i,j) = 0; % set val that colse to 0 to 0
            if (magnitude(i, j-1) > magnitude(i, j) || magnitude(i, j+1) > magnitude(i, j))
                im_rec(i, j) = 0;
            end
        elseif (angl(i, j) >= 22.5) || (angl(i, j) < 67.5) % 45 degree
            angl(i,j) = 45;
            if (magnitude(i-1, j-1) > magnitude(i, j) || magnitude(i+1, j+1) > magnitude(i, j))
                im_rec(i, j) = 0;
            end
        elseif (angl(i, j) >= 67.5) || (angl(i, j) < 112.5) % 90 degree
            angl(i,j) = 90;
            if (magnitude(i-1, j-1) > magnitude(i, j) || magnitude(i+1, j+1) > magnitude(i, j))
                im_rec(i, j) = 0;
            end
        elseif (angl(i, j) >= 112.5) || (angl(i, j) < 157.5) % 45 degree
            angl(i,j) = 135;
            if (magnitude(i-1, j-1) > magnitude(i, j) || magnitude(i+1, j+1) > magnitude(i, j))
                im_rec(i, j) = 0;
            end
        else
        end
    end
img1 = im_rec;
end
end

    
                
        
        
