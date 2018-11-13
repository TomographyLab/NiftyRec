function image_out = et_resize(image_in, nx, ny, nz)

    [sz,sx,sy] = size(image_in);
    [x,y,z]    = meshgrid(1:1:sx,1:1:sy,1:1:sz);
    [xi,yi,zi] = meshgrid(1:sx/nx:sx,1:sy/ny:sy,1:sz/nz:sz);
    image_out  = interp3(x,y,z,image_in,xi,yi,zi);
end