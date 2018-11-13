function backprojection = et_backproject_irt(sinogram, cameras, attenuation, psf, use_gpu, background, background_attenuation, truncate_negative_values)

%ET_BACKPROJECT_IRT
%    Backprojection function for Emission Tomographic reconstruction
%    Same result as 'et_backproject', but uses IRT 
%
%Description
%    Function for backprojection of activity from detector space.
%
%    IMAGE = ET_BACKPROJECT_IRT(SINOGRAM, CAMERAS, ATTENUATION, PSF)
%
%    SINOGRAM is a 2D or 3D sinogram.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    ATTENUATION specifies attenuation coefficients. 
%
%    PSF is a Depth-Dependent Point Spread Function.
%
%Reference
%    IRT reconstruction toolbox, Fessler. 
%
%Example
%    N = 128;
%    n_cameras = 100;
%    use_gpu = 1;
%    background = 0;
%    background_attenuation=0;
%    sinogram = ones(N,N,n_cameras);
%    attenuation = zeros(N,N,N);
%    PSF = ones(11,11,N);
%    cameras = [0:pi/n_cameras:pi-pi/n_cameras]';
%    image = et_backproject(sinogram,cameras,attenuation,PSF,use_gpu,background,background_attenuation);
%
%See also
%    ET_BACKPROJECT, ET_PROJECT_IRT, ET_MLEM_DEMO, ET_OSEM_DEMO
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK



if not(exist('attenuation'))
    attenuation = 0;
end

N_x = size(sinogram,1);
N_y = size(sinogram,1);
N_z = size(sinogram,2);
N_projections = size(sinogram,3); 
if size(cameras,2)==3
    rotation_angle = 180/pi*(((cameras(N_projections,2)-cameras(1,2))/N_projections)*(N_projections+1));
    initial_rotation_angle = 180/pi*cameras(1,2);
else
    rotation_angle = 180/pi*(((cameras(N_projections)-cameras(1))/N_projections)*(N_projections+1));
    initial_rotation_angle = 180/pi*cameras(1);
end

f.chat = 0;
f.option = {};
f.test = '3s';
dx=1;
ig = image_geom('nx', N_x, 'ny', N_y, 'nz', N_z, 'dx', dx, 'dz', dx);
if isscalar(psf)
    f.psfs = '-';
else
    f.psfs = '/tmp/t,psfs.fld';
    fld_write(f.psfs, psf)
end
f.blur = ',fft'; 
f.na = N_projections;
f.dir = test_dir;
if attenuation==0
    f.mumap = '-';
else
    f.mumap = [f.dir 'mumap.fld'];
    fld_write(f.mumap, attenuation); 
end
f.sfilter = 1; 
f.sys_type = '3s@%g,%g,%g,%g,%g,%d%s@%s@%s@-%d,%d,%d';
f.sys_type = sprintf(f.sys_type, ig.dx, ig.dx, ig.dz, rotation_angle,initial_rotation_angle, f.sfilter, f.blur, f.mumap, f.psfs, ig.nx, ig.nz, f.na); 
f.option = {f.option{:}, 'chat', f.chat, 'checkmask', 0};
G = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, f.option{:}, 'nthread', jf('ncore'));
backprojection2 = G'*sinogram;
backprojection = zeros(N_x,N_z,N_y);
for i=1:N_z
    backprojection(:,i,:)=backprojection2(:,:,i);
end





