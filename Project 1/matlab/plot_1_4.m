x1 = 10;
x2 = 20;
x3 = 30;
theta = 0:pi/100:2*pi;
rho1 = x1*cos(theta) + x1*sin(theta);
rho2 = x2*cos(theta) + x2*sin(theta);
rho3 = x3*cos(theta) + x3*sin(theta);
plot(theta, rho1, theta, rho2, theta, rho3)

axis([0 5.0 -50 50])