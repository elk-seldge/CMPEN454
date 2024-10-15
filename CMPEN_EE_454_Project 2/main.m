% 3D-TO-2D & 2D-TO-3D & triangulation f*X/Z f*Y/Z
load('Subject4-Session3-Take4_mocapJoints.mat', 'mocapJoints');
load('vue2CalibInfo.mat', 'vue2');
load('vue4CalibInfo.mat', 'vue4');

% 2.1 
mocapFnum = 1000;  %mocap frame number 1000
x = mocapJoints(mocapFnum,:,1);     %array of 12 X coordinates 
y = mocapJoints(mocapFnum,:,2);     %  Y coordinates
z = mocapJoints(mocapFnum,:,3);     %  Z coordinates 
conf = mocapJoints(mocapFnum,:,4);   %confidence values

% 2.2 
% focal length
v2_fcl = vue2.foclen;
v4_fcl = vue4.foclen;
% principle point
v2_cx = vue2.prinpoint(1);
v2_cy = vue2.prinpoint(2);
v4_cx = vue4.prinpoint(1);
v4_cy = vue4.prinpoint(2);
% 
v2_Sx = 1;
v2_Sy = 1;
v4_Sx = 1;
v4_Sy = 1;
% rotation matrix
v2_R = vue2.Rmat;
v4_R = vue4.Rmat;
% offset
v2_Ox = vue2.position(1);
v2_Oy = vue2.position(2);
v2_Oz = vue2.position(3);

v4_Ox = vue4.position(1);
v4_Oy = vue4.position(2);
v4_Oz = vue4.position(3);
% verify
flm_to_pxl = [ 1/v2_Sx   0        v2_cx ;
                0        1/v2_Sy  v2_cy ;
                0        0        1     ];
    
p_p = [ v2_fcl  0       0   0 ;
        0       v2_fcl  0   0 ;
        0       0       1   0 ];

disp("(flm2pxl) * (perspective projection):");
disp(flm_to_pxl * p_p);
disp("K mat:");
disp(vue2.Kmat);
% disp("K mat == flm2pxl * pp");
% rec = (flm_to_pxl * p_p == vue2.Kmat);
% disp(rec);

rttn = [ v2_R(1,1) v2_R(1,2) v2_R(1,3)  0;
         v2_R(2,1) v2_R(2,2) v2_R(2,3)  0;
         v2_R(3,1) v2_R(3,2) v2_R(3,3)  0;
         0          0        0          1];

offst = [1  0   0   -v2_Ox;
         0  1   0   -v2_Oy;
         0  0   1   -v2_Oz;
         0  0   0       1];

disp("rotation * offset:")
disp(rttn * offst);
disp("P mat:");
disp(vue2.Pmat);

% 2.3
%initialization of VideoReader for the vue video. 
%YOU ONLY NEED TO DO THIS ONCE AT THE BEGINNING
filenamevue2mp4 = 'Subject4-Session3-24form-Full-Take4-Vue2.mp4'; 
filenamevue4mp4 = "Subject4-Session3-24form-Full-Take4-Vue4.mp4";
vue2video = VideoReader(filenamevue2mp4);
vue4video = VideoReader(filenamevue4mp4);
%now we can read in the video for any mocap frame mocapFnum.
%the (50/100) factor is here to account for the difference in frame 
%rates between video (50 fps) and mocap (100 fps).
vue2video.CurrentTime = (mocapFnum-1)*(50/100)/vue2video.FrameRate; 
vue4video.CurrentTime = (mocapFnum - 1)*(50/100)/vue4video.FrameRate;
vid2Frame = readFrame(vue2video);
vid4Frame = readFrame(vue4video);

figure(1); % create image window
image(vid2Frame); % arr to image, display
hold on; % do not flash, keep axis

% 2d-vue2
v2_2d = zeros(3,12);
for joint_num = 1:12
    [vue2_u, vue2_v] = point3D_to_pixel2D(v2_Sx, v2_Sy, v2_cx, v2_cy, v2_R, v2_Ox, v2_Oy, v2_Oz, x(joint_num), y(joint_num), z(joint_num));
    v2_2d(1,joint_number) = vue2_u;
    v2_2d(2,joint_number) = vue2_v;
    v2_2d(3,joint_number) = 1;
    plot(vue2_u, vue2_v, 'r*', 'MarkerSize', 5, 'LineWidth', 2);
end

xlim = get(gca,'XLim'); % gca ---> Cartesian coordinate system

m = (1.1266667 - v2_2d(2,1)) / ((-147.893) - v2_2d(1,1));
n = 1.1266667 - (-147.893) * m;

y1 = m * xlim(1) + n;
y2 = m * xlim(2) + n;
line([xlim(1) xlim(2)], [y1 y2]);

figure(1); % create image window
image(vid2Frame); % arr to image, display
hold on; % do not flash, keep axis

% 2d-vue4
v4_2d = zeros(3,12);
for joint_num = 1:12
    [vue4_u, vue4_v] = point3D_to_pixel2D(v4_Sx, v4_Sy, v4_cx, v4_cy, v4_R, v4_Ox, v4_Oy, v4_Oz, x(joint_num), y(joint_num), z(joint_num));
    v4_2d(1,joint_number) = vue4_u;
    v4_2d(2,joint_number) = vue4_v;
    v4_2d(3,joint_number) = 1;
    plot(vue4_u, vue4_v, 'r*', 'MarkerSize', 5, 'LineWidth', 2);
end

% 2.4 Triangulation
R_r = v2_R;
R_l = v4_R;

C_r = [v2_Ox;
       v2_Oy; 
       v2_Oz];
C_l = [v4_Ox;
       v4_Oy;
       v4_Oz];
% internal
K_r = vue2.Kmat;
K_l = vue4.Kmat;

v2_3d = zeros(3,12);
% 2dto3d-vue2
for joint_number = 1:12
    P_r = [v2_2d(1,joint_number)
           v2_2d(2,joint_number)
           1];
    P_l = [v4_2d(1,joint_number)
           v4_2d(2,joint_number)
           1];
    [u,v,w] = pixel2D_to_point3D(R_r, R_l, C_r, C_l, K_r, K_l, P_r, P_l);
    v2_3d(1,joint_number) = u;
    v2_3d(2,joint_number) = v;
    v2_3d(3,joint_number) = w;
end
% Display calculated 3d point for 1st joint
disp("Inverse Projection: ");
disp(v2_3d(1,1) + ", " +  v2_3d(2,1) + ", " +  v2_3d(3,1));

% Display actual 3d point for 1st joint
disp("Original Point: ");
disp(x(1) + ", " + y(1) + ", " + z(1));

% 2.5 
euc_L2 = 0;
err = 0;
euc_sum = 0;

for joint_number = 1:12
    % compare points in v2_3d to original x, y, and z
    dx = x(joint_number) - v2_3d(1,joint_number);
    dy = y(joint_number) - v2_3d(2,joint_number);
    dz = z(joint_number) - v2_3d(3,joint_number);
    
    % calculated for sum of	L2 distances
    euc_L2 = sqrt(dx^2 + dy^2 + dz^2);
    euc_sum = euc_sum + euc_L2;
end

err = euc_sum / 12;
disp("error:");
disp(err);

% 2.6 epippolar
R = (vue2.Rmat)^(-1) * (vue4.Rmat);
disp(rank(R));
S = [0                  -(v4_Oz - v2_Oz)    (v4_Oy - v2_Oy);
     (v4_Oz - v2_Oz)    0                   -(v4_Ox - v2_Ox);
     -(v4_Oy - v2_Oy)   (v4_Ox - v2_Ox)                   0];

E = R * S;
disp(rank(E));

e_linear = [0 0 0];
e_v2 = E \ e_linear.';
e_v4 = -(linsolve(E,e_linear.'));
disp("rotation matrix:");
disp(R);
disp("essential matrix:");
disp(E);
disp("v2 epipole:");
disp(e_v2);
disp("v4 epipole:");
disp(e_v4);

% 3D to 2D func
function [u,v] = point3D_to_pixel2D(Sx, Sy, Cx, Cy, fcl, R_arr, Ox, Oy, Oz, Pu, Pv, Pw)

    flm_to_pxl = [ 1/Sx 0 Cx ;
                   0 1/Sy Cy ;
                   0 0 1 ];
    
    p_p = [ fcl 0 0 0;
            0 fcl 0 0;
            0 0 1 0];
                          
    rotation = [ R_arr(1,1) R_arr(1,2) R_arr(1,3) 0;
                 R_arr(2,1) R_arr(2,2) R_arr(2,3) 0;
                 R_arr(3,1) R_arr(3,2) R_arr(3,3) 0;
                 0          0          0          1];

    offset = [1 0 0 -Ox;
              0 1 0 -Oy;
              0 0 1 -Oz;
              0 0 0   1];
          
    world_coord = [Pu;
                   Pv;
                   Pw;
                    1];
                  
    pxl_locat = flm_to_pxl * p_p * rotation * offset * world_coord;
        
    u = pxl_locat(1) / pxl_locat(3); 
    v = pxl_locat(2) / pxl_locat(3); 
end

% 2D to 3D func
function [u,v,w] = pixel2D_to_point3D(R_r, R_l, C_r, C_l, K_r, K_l, P_r, P_l)
    % calculate direction vectors
    U_r = (R_r.')*(K_r^(-1))*P_r;
    U_l = (R_l.')*(K_l^(-1))*P_l;
    
    % calculate normalized direction vectors
    U_r_hat = U_r / norm(U_r);
    U_l_hat = U_l / norm(U_l);
    
    % calculate cross product norm
    U_x_hat = cross(U_l_hat, U_r_hat) / norm(U_r_hat + U_l_hat);
    
    % solve linear system
    syms a b d
    matrix_eqn = a*U_l_hat + d*U_x_hat - b*U_r_hat == C_r - C_l;
    
    % condense equation into an augmented matrix:  X | y
    [X, y] = equationsToMatrix(matrix_eqn, [a b d]);

    % solve augmented matrix
    Z = linsolve(X, y);
    a = Z(1);
    b = Z(2);
    % d = Z(3);
    
    % calculate right 3d point
    P_3d_r = C_r + b*U_r_hat;
    
    % calculate left 3d point
    P_3d_l = C_l + a*U_l_hat;
    
    % calculate midpoint of right / left 3d point
    P_3d = (P_3d_r + P_3d_l) / 2;
    
    % cast P_3d to double
    P_3d = double(P_3d);

    % output solution
    u = P_3d(1);
    v = P_3d(2);
    w = P_3d(3);

end
