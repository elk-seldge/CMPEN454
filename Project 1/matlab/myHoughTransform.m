function [H, rhoScale, thetaScale] = myHoughTransform(Im, threshold, rhoRes, thetaRes)
%Your implementation here
% rhoRes = 2; thetaRes  = pi/90;threshold = 0.03; 
% rho need to be large enough
% accomodate all lines that could lie in an image
% rhoMax = 1 * sqrt(row * row + col * col)
%thetaMax = 360 degree

% voteMatrix = zeros(colRho, colThres)
% for 
% if pixel > threshold
% find corresponding rho and theta
% voteMatrix[i, j] += 1

% 1.  Record vote for each possible line on which each edge point lies.
% 2.  Look for lines that get many votes.
% rhoScale and thetaScale are arrays of rho val and theta val to generate H
% using Hough transform

[rows, cols] = size(Im);
rhoMax =  1 * sqrt(rows * rows + cols * cols); % 800
rhoScale = -rhoMax:rhoRes:rhoMax; % rhoRes = 2 rhoScale = [1, 801]
thetaScale = -pi/2: thetaRes * thetaRes: pi/2; % thetaRes resolution on theta axis pi/90 thetaScale =[1, 5157]
H = zeros(size(rhoScale,2), size(thetaScale, 2)); % numel return the number of elements
realThetaScale = (90/pi)*thetaRes;
for x = 1: rows
    for y = 1: cols
        if (Im(x, y) > threshold)
            for theta = thetaScale
                    tht = (theta/pi)*180;
                    rho = y * cosd(tht)+ x * sind(tht);
                    H(floor((rho+rhoMax)/rhoRes)+1,floor((tht+90)/realThetaScale)+1) = H(floor((rho+rhoMax)/rhoRes)+1,floor((tht+90)/realThetaScale)+1) + 1;
            end
        end
    end
end

end
        
        