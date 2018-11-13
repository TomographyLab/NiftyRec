
%% Test irotate

%% Parameters
N             = 128;
N_rotations   = 120;
N_spheres     = 40; 

GPU           = 1;

%% Make phantom
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
phantom = zeros(N,N,N); 
for i=1:N_spheres
    radius = N/16 + N/8*rand(); 
    foreground = 0.2+rand();
    background = 0.01;
    x = randint(1,1,N); y = randint(1,1,N); z = randint(1,1,N); 
    phantom = phantom + et_spherical_phantom(N,N,N, radius, foreground, background, x, y, z) .* mask;
end

%% Rotate along random trajectory and irotate
dx = (0.1+0.9*rand())/180*pi; % 0.1 to 1 degrees 
dy = (0.1+0.9*rand())/180*pi; 
dz = (0.1+0.9*rand())/180*pi; 

figure; 
for i=1:N_rotations
    rotation = et_rotate(phantom,i*[dx,dy,dz],[(N+1)/2,(N+1)/2,(N+1)/2],GPU); 
    irotation = et_irotate(rotation,i*[dx,dy,dz],[(N+1)/2,(N+1)/2,(N+1)/2],GPU); 
    subplot(1,3,1); image(N_spheres*phantom(:,:,64)); axis square; title('phantom')
    subplot(1,3,2); image(N_spheres*rotation(:,:,64)); axis square; title('rotation')
    subplot(1,3,3); image(N_spheres*irotation(:,:,64)); axis square; title('i-rotation')
    pause(0.05); 
end


