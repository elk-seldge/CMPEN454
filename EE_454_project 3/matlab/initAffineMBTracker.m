function [affineMBContext] = initAffineMBTracker(img, rect)
% img is a greyscale image with a bouunding box rect
% affineMBContext is a Matlab structure that contains the Jacobian of the
% affine warp with respect to the 6 affine warp parameters and the inverse
% of the approximated Hessian matrix (J and H^-1)

x=rect(1);
y=rect(2);
w=rect(3);
h=rect(4);

T = double(img(y:y+h-1, x:x+w-1))/255.0;
[dTx,dTy] = gradient(T);

% flatten to column vectors
fprintf('size T');
disp(size(T));
T = T(:);
dT = [dTx(:), dTy(:)];

[X,Y] = meshgrid(x:x+w-1, y:y+h-1);
fprintf('size X');
disp(size(X));
X=X(:);
Y=Y(:);

%[Y,X] = meshgrid(y:y+h-1,x:x+w-1);
%    X=X(:);
%    Y=Y(:);

o=ones(size(X));

J = [dT dT dT].*[X X Y Y o o];

affineMBContext.J = J;
affineMBContext.X = X;
affineMBContext.invH = (J'*J)^(-1);

end

