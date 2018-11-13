function RES = et_spherical_phantom(x_size, y_size, z_size, radius, in_value, out_value, x_center, y_center, z_center)

%ET_SPHERICAL_PHANTOM
%    Creates a sphere of uniform activity in a uniform background
%
%Description
%    Function that creates a 3D activity image with uniform activity in a spherical region
%    and uniform background. 
%
%    IMAGE = ET_SPHERICAL_PHANTOM(X_S, Y_S, Z_S, RADIUS, IN_VALUE, OUT_VALUE, X_C, Y_X, Z_C)
%
%    X_S, Y_S, Z_S specify the size of the image.
%
%    RADIUS specifies the radius of the sphere.
%
%    IN_VALUE specifies the activity within the sphere
%
%    OUT_VALUE specifies the activity outside the sphere
%
%    X_C, Y_X, Z_C specify the center of the sphere in pixels
%
%Example
%    N = 128;
%    activity = 100;
%    background = 10;
%    phantom = et_spherical_phantom(N,N,N,N/8,activity,background,N/2,N/2,N/2);
%
%See also
%   ET_BACKPROJECT, ET_BACKPROJECT, ET_MAPEM_STEP
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


if not(exist('x_size','var'))
    x_size = 64;
end
if not(exist('y_size','var'))
    y_size = 64;
end
if not(exist('z_size','var'))
    z_size = 64;
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
if not(exist('radius','var'))
    radius = min([x_size,y_size,z_size])/8;
end
if not(exist('out_value','var'))
    out_value = 0;
end
if not(exist('in_value','var'))
    in_value = 1;
end

RES = zeros(x_size,y_size,z_size)+out_value;

for x = 1 : x_size
    for y = 1 : y_size
        for z = 1 : z_size
            d = sqrt( (x-x_center)^2 + (y-y_center)^2  + (z-z_center)^2 ) ;
            if ( d < radius )
                RES(x,y,z) = in_value ;
            end
        end
    end
end

