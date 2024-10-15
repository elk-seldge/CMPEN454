function [u, v] = LucasKanade(It, It1, rect)

[m0,n0] = size(It);
%[m1,n1] = size(It1);

x = rect(1);
y = rect(2);
w = rect(3);
h = rect(4);
t = y;
l = x; 
b = y+h-1; 
r = x+w-1;

dp = ones(6,1) * 50; % * num does influence the result. 
thrshld = 3;
itertime = 30;
p = zeros(6,1);
for i=1:itertime
    if norm(dp) < thrshld
        break;
    end
    Wx = WP(p);
    % 1. I(W(x;p))
    warpedIt1 = warpH(It1, Wx, [m0, n0], 0);

    % 2. T(x)-I(W(x;p))
    I = double(warpedIt1(t:b,l:r));
    T = double(It(t:b,l:r));
    E = (T-I);
    
    % 3. ∇I(x)
    [d_x, d_y] = gradient(I); % ∇I(x)
    dx = d_x(:); % m by 1
    dy = d_y(:); % m by 1

    % 4. dW/dp and compute A
    [xs, ys] = meshgrid(t:b,l:r);
    xs = xs(:); % m by 1
    ys = ys(:); % m by 1
    
    J = [xs.*dx, xs.*dy, ys.*dx, ys.*dy, dx, dy]; % m by 6
    
    % 5. H --> A^T*A
    H = J'*J; % 6 by 6
    
    % 6. dp main problem may locate in
    dp = H\(J'*E(:)); % E --> b
    % 7. p->p+dp
    p = p - dp;
    disp(p)
end

u = p(5);
v = p(6);
fprintf('u: %d, v: %d\n',u, v);
end 