function volume = et_spiral_phantom(x_size, y_size, z_size, foreground_value, background_value, n_spheres, min_radius_sphere, max_radius_sphere, radius_spiral, axis, x_center, y_center, z_center)

%ET_SPIRAL_PHANTOM
%    Spiral phantoms. 
%
%Description
%    This function creates a volumetric image with uniform background and 
%    and uniform value in spherical regions displaced on a circle. 
%
%    IMAGE =
%    ET_SPIRAL_PHANTOM(X_SIZE,Y_SIZE,Z_SIZE,FOREGROUND_VALUE,BACKGROUND_VALUE,N_SPHERES,MIN_RADIUS_SPHERE,MAX_RADIUS_SPHERE,RADIUS_SPIRAL,AXIS,CENTER_X,CENTER_Y,CENTER_Z)
%
%    X_SIZE, Y_SIZE, Z_SIZE specify the size of the image (default 128,32,128)
%
%    FOREGROUND_VALUE specifies the value of the foreground voxels (default 1)
%
%    BACKGROUND_VALUE specifies the value of the background voxels (default 0)
%
%    N_SPHERES specifies the number of spheres (default 8)
%
%    MIN_RADIUS_SPHERE radius of the smallest sphere (default 1)
%
%    MAX_RADIUS_SPHERE radius of the largest sphere (default 15)
%
%    RADIUS_SPIRAL specifies the radius of the spiral in units of voxels (default 45)
%
%    AXIS specifies the axis of rotation: 1->x 2->y 3->z (default 2)
%    CENTER_X, CENTER_Y,CENTER_Z specify the center of the sphere in units
%    of voxels (default X_SIZE/2,Y_SIZE/2,Z_SIZE/2)
%
%
%Example
%    N = 128;
%    n_spheres = 8;
%    foreground = 100;
%    background = 10;
%    phantom = et_spiral_phantom(N,N,N,foreground,background,n_spheres,2,10,40,2); 
%    projection = et_project(phantom,linspace(0,pi,180)',0,0,1);
%
%See also
%   ET_SPHERICAL_PHANTOM
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


if not(exist('x_size','var'))
    x_size = 128;
end
if not(exist('y_size','var'))
    y_size = 32;
end
if not(exist('z_size','var'))
    z_size = 128;
end
if not(exist('foreground_value','var'))
    foreground_value = 1;
end
if not(exist('background_value','var'))
    background_value = 0;
end
if not(exist('n_spheres','var'))
    n_spheres = 8;
end
if not(exist('min_radius_sphere','var'))
    min_radius_sphere = 1;
end
if not(exist('max_radius_sphere','var'))
    max_radius_sphere = 15;
end
if not(exist('radius_spiral','var'))
    radius_spiral = 45;
end
if not(exist('axis','var'))
    axis = 2;
end
if not(exist('x_center','var'))
    x_center = round(x_size/2);
end
if not(exist('y_center','var'))
    y_center = round(y_size/2);
end
if not(exist('z_center','var'))
    z_center = round(z_size/2);
end


volume = zeros(x_size,y_size,z_size)+background_value; 
i=0;
for radius_sphere = linspace(min_radius_sphere,max_radius_sphere,n_spheres)
    theta = 2*pi / n_spheres * i;
    loc_sphere = [cos(theta),sin(theta);-sin(theta),cos(theta)]*radius_spiral*[0.5*sqrt(2),0.5*sqrt(2)]'; 
    if axis == 1
        x_sphere = x_center; 
        y_sphere = loc_sphere(1)+y_center;
        z_sphere = loc_sphere(2)+z_center; 
    elseif axis == 3
        x_sphere = loc_sphere(1)+x_center; 
        y_sphere = loc_sphere(2)+y_center;
        z_sphere = z_center;         
    else
        x_sphere = loc_sphere(1)+x_center; 
        y_sphere = y_center;
        z_sphere = loc_sphere(2)+z_center;         
    end
    volume = volume + et_spherical_phantom(x_size,y_size,z_size,radius_sphere,foreground_value,0,x_sphere,y_sphere,z_sphere); 
    i=i+1;
end


