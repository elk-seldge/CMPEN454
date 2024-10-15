function [W_out] = affineMBTracker(img, tmp, rect, W_in, context)
% img is a greyscale image
% tmp is the template image
% rect is the bounding box
% W_in is the previous frame
% context is the precomputed J and H^-1 matrices
%
% W_out should be 3 by 3 that contains the new affine warp matrix updated so
% that it aligns the CURRENT frame with the TEMPLATE

itertime = 30;
x = rect(1);
y = rect(2);
w = rect(3);
h = rect(4);

[X,Y] = meshgrid(x:x+w-1, y:y+h-1);
X = X(:);
Y = Y(:);
o = ones(size(X));

T = tmp(:);
invH = context.invH;
J = context.J;
W_out = W_in;

for i = 1: itertime
    C = W_out * [X';
                 Y';
                 o'];
    C = C ./ C(3,3);
    I = img(sub2ind(size(img), floor(C(2,:)'), floor(C(1,:)')));
    I = I .* (sum(T(:)) / sum(I(:)));
    D = (double(I)-double(T))/255.0;
    dp = invH * J' * D;
    W_out = W_out / WP(dp);
    diff = sum(abs(D))/length(D);
end
% fprintf("iter %i, dp %i, diff %i,",i, norm(dp), diff)
end
