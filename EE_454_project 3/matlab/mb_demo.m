clear
clc

x_center = 240;
y_center = 192;
width = 210;
height = 170;
l = x_center-floor(width/2); %
t = y_center-floor(height/2); %
tracker = [l t width height]; %

video = VideoWriter('../results/car_mb_2.avi'); % result
open(video); %

% Initialize the tracker
figure;


% TODO run the Matthew-Baker alignment in both landing and car sequences
prev_frame = imread('../data/car/frame0020.jpg'); % 
template = prev_frame(t:(t+height-1),l:(l+width-1));   % TODO
imshow(prev_frame);
%implement your code here
%---------------------------------------

%Win = WP([0 0 0 0 0 0]);
Win = [1 0 0 ; 0 1 0 ; 0 0 1 ] ;
context = initAffineMBTracker(prev_frame, tracker);
%---------------------------------------
 

% Start tracking
new_tracker = tracker;
for i = 21:280
    imgdir = sprintf('../data/car/frame%04d.jpg', i);
    if (~exist(imgdir,'file'))
        continue;
    end
    im = imread(imgdir);
    Wout = affineMBTracker(im, template, tracker, Win, context);

    
    %implement your code here
    %---------------------------------------
    new_c = Wout * [tracker(1) ; tracker(2) ; 1];
%     new_c = Wout * [tracker(1) + tracker(3)/2;
%                     tracker(2) + tracker(4)/2;
%                     1];
%    new_tracker = [new_c(1)/new_c(3) - width/2, new_c(2)/new_c(3) - height/2, width, height];
    new_tracker = [new_c(1) new_c(2) tracker(3) tracker(4)];
%---------------------------------------

    clf;
    hold on;
    imshow(im);   
    
    disp(new_tracker);
    rectangle('Position', new_tracker, 'EdgeColor', [1 1 0]);
    title(num2str(i));
    drawnow;

    F = getframe;
    writeVideo(video,F);

end

close(video);
