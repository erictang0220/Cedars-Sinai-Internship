function [] = Coil_Brain_Plot_9var(input, mask, bodyMask)
%2016/12/28 Randy Yang
% Bz: Z-direction magnetic field 
%-----------r(meter),theta(degree)_rotate againsty,phi(degree)_rotate against z: coil center location to the center of the hemsphere-
%-----------a(meter): radius of the coil               
%-----------abratio: elliptical ratio
%-----------angX(degree):  angle against x of the tangent plane
%-----------angY(degree):  angle against y the tangent plane
%-----------angZ(degree):  angle against z of the tangent plane
%% parameter setup

% physical phantom
% FOVx = 0.400; % FOV of the B0 maps to shim(m) in x, y direction
% FOVy = 0.300;

% original
FOVx = 0.380; 
FOVy = 0.380;

% digital phantom
% FOVx = 0.45; 
% FOVy = 0.45;

% physical phantom
% nx0 = 112; % number of slices (array dimension) in x, y, z direction
% ny0 = 84;
% nz0 = 36;

%original
nx0 = 192; 
ny0 = 192; 
nz0 = 96;

% digital phantom
% nx0 = 512; 
% ny0 = 512;
% nz0 = 40;

% digital phantom 3 only
% nz0 = 50;

resX = FOVx/nx0; % resolution (voxel x, y dimension) of the B0 maps
resY = FOVy/ny0;

% physical phantom
% height = 0.2361; % FOV in z-direction

% original
height = 0.192; 

% digital phantom
% height = 0.2; 

% digital phantom 3 only
% height = 0.25;

th = height/nz0; % thickness, aka z-resolution

%I = 1; % reference current in each coil (A)

% physical coordinates
% x,y,z will have the same dimension
% use FOV directly for x and y, no extra padding
[x,y,z] = meshgrid(-(FOVy-resY)/2:resY:(FOVy-resY)/2,...
    -(FOVx-resX)/2:resX:(FOVx-resX)/2, -(height-th)/2:th:(height-th)/2);

% Get the middle portion of the Mask
% mask = mask(end/4+1:end*3/4,:,:);
% bodyMask = bodyMask(end/4+1:end*3/4,:,:);

%% correct for an offset of the center
    
% [nx ny nz]= size(x);
% nxi = round(nx/2)-round(nx0/2); nxf = round(nx/2)+round(nx0/2)-1;
% nyi = round(ny/2)-round(ny0/2); nyf = round(ny/2)+round(ny0/2)-1;
% nzi = round(nz/2)-round(nz0/2); nzf = round(nz/2)+round(nz0/2)-1;

% offset the prostate to the origin
% offset for z (in meters)

% original and physical phantom
rhz = -0.06;

% digital phantom
% rhz = -0.015;
% rhz = -0.025;

% offset for x (in meters)

% original
rhx = 0.005;

% physical phantom
% rhx = 0;

% digital phantom
% rhx = 0.05;
% rhx = 0.025;
% rhy = 0.025;

z = z + rhz;
x = x + rhx;
% y = y + rhy;

% x = x(nxi-offsetx:nxf-offsetx,nyi:nyf,nzi-offsetz:nzf-offsetz);
% y = y(nxi-offsetx:nxf-offsetx,nyi-offsety:nyf-offsety,nzi-offsetz:nzf-offsetz);
% z = z(nxi-offsetx:nxf-offsetx,nyi-offsety:nyf-offsety,nzi-offsetz:nzf-offsetz);

Bz = zeros([size(x,1) size(x,2) size(x,3) size(input,1)]);

Bzc = zeros([size(x,1) size(x,2) size(x,3)]);
coils = zeros([size(x,1) size(x,2) size(x,3)]);

% Mask the prostate
brainx = x(mask);
brainy = y(mask);
brainz = z(mask);

% Mask the body
bodyx = x(bodyMask);
bodyy = y(bodyMask);
bodyz = z(bodyMask);

figure; scatter3(brainx(1:10:end), brainy(1:10:end), brainz(1:10:end), 'r');
hold on;
s = scatter3(bodyx(1:10:end), bodyy(1:10:end), bodyz(1:10:end), 'b');
alpha(s, 0.1);

%% Ellipse field generation
N=size(input,1);
for coil=1:N
    % the center coordinate of the circle 
    xc=input(coil,1); 
    yc=input(coil,2); 
    zc=input(coil,3); 
    
    % radius (meter)
    a=input(coil,4); 
    % b = abratio * a
    abratio=input(coil,5); 
    
    % angle of rotation (degree) around x,y,z axis
    angX=input(coil,6); 
    angY=input(coil,7);
    angZ=input(coil,8);
    
    I=input(coil,9);
    
    % from degree to radian
    angX=angX/180*pi;
    angY=angY/180*pi;
    angZ=angZ/180*pi;

%------------------draw the circle on x-y plane-------------------
    b = a*abratio;

    x1 = a*cos((0:359)*pi/180);
    x2 = a*cos((1:360)*pi/180);

    z1 = zeros(1,size(x1,2));
    z2 = z1;

    y1 = b*sin((0:359)*pi/180);
    y2 = b*sin((1:360)*pi/180);

    coordi1 = [x1; y1; z1];
    coordi2 = [x2; y2; z2];
    %-------------------------


    %% 1 raws
    
    % for rot3, rot2, rot1, rotation is positive when CCW viewed from the
    % axis
    
    
    % rot1: rotation about the x-axis
    rot1 = [1 0 0;
            0 cos(angX) -sin(angX); 
            0 sin(angX) cos(angX)];
    % rot2: rotation about the y-axis
    rot2 = [cos(angY) 0 sin(angY); 
            0 1 0;  
            -sin(angY) 0 cos(angY)];
    % rot3: rotation about the z-axis
    rot3 = [cos(angZ) -sin(angZ) 0; 
            sin(angZ) cos(angZ) 0; 
            0 0 1];
        
    % rotation matrix: zyx!!! Order matters!!!
    coordi1_2 = rot3*rot2*rot1*coordi1;
    % coordi2_2 = rot1*rot2*rot3*rot4*rot5*coordi2+coffset;

    coordi1_1 =  coordi1_2;

    % x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    % y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    % z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
        
    % move to the center in the end
    x1f = coordi1_2(1,:) + xc;
    y1f = coordi1_2(2,:) + yc;
    z1f = coordi1_2(3,:) + zc;
    
    scatter3(x1f,y1f,z1f)


end
    %seperate Bz for each coil for parrallization

    
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])



end



%Bz_f = Bzsum(nxi-offsety:nxf-offsety,nyi:nyf,nzi-offsetz:nzf-offsetz);
