
% TT_DEMO_MLEM_PARALLEL_LARGE
%     NiftyRec Demo: iterative reconstruction transmission tomography with
%     parallel rays - large volume
%
%See also
%   TT_DEMO_PROJECT_PARALLEL_LARGE
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK

% The attenuation is loaded from file 'attenuation_phantom_file_nifti'; 
% The atteunation is then cropped and padded in order to fill the volume 
% of size Nx,Ny,Nz specified below. 
% The only constraint on Nx,Ny,Nz is that the largest of Nx,Ny,Nz must be larger 
% than the largest volume size (an error is returned if this condition is
% not verified). The intended use is: Nx=Nz='multiple of 64 just about larger than the
% largest dimension of the volume loaded from file' and 'Ny=desired number
% of slides in axial direction'. 
% Volume is initially optionally rotated around one of its axis by 90 deg 
% in order to select the desired axial direction. 

% Load phantom attenuation data, pad, rotate, crop
attenuation_phantom_file_nifti = 'attenuation_01.nii'; 
Nx = 448; 
Ny = 64; 
Nz = 448; 
rotate_axis = 3;   % 0 -> no rotation; 1 -> x; 2 -> y; 3 -> z

fprintf('Loading phantom ...\n');
N = max([Nx,Ny,Nz]);
attenuation = et_load_nifti(attenuation_phantom_file_nifti);
attenuation = single(attenuation.img)*5e-7;
attenuation(attenuation<=0)=0;
[nx,ny,nz] = size(attenuation);
fprintf('Padding ...\n');

attenuation_phantom = zeros(N,N,N); 
attenuation_phantom(floor((N-nx)/2)+1:floor((N-nx)/2)+nx,floor((N-ny)/2)+1:floor((N-ny)/2)+ny,floor((N-nz)/2)+1:floor((N-nz)/2)+nz) = attenuation;
clear attenuation

fprintf('Rotating ...\n'); 
attenuation_phantom2 = attenuation_phantom; 
if rotate_axis ~= 0
    for i =1:Nz
        if rotate_axis == 1
            attenuation_phantom(:,i,:) = attenuation_phantom2(:,:,i);
        elseif rotate_axis == 2
            attenuation_phantom(i,:,:) = attenuation_phantom2(:,i,:);            
        elseif rotate_axis == 3            
            attenuation_phantom(:,:,i) = attenuation_phantom2(i,:,:);            
        end
    end
end
clear attenuation_phantom2

fprintf('Cropping ...\n');
if Nx < N
    attenuation_phantom = attenuation_phantom(floor((N-Nx)/2)+1:floor((N-Nx)/2)+Nx,:,:); 
end
if Ny < N
    attenuation_phantom = attenuation_phantom(:,floor((N-Ny)/2)+1:floor((N-Ny)/2)+Ny,:); 
end
if Nz < N
    attenuation_phantom = attenuation_phantom(:,:,floor((N-Nz)/2)+1:floor((N-Nz)/2)+Nz); 
end

% Set the parameters of the simulation and reconstruction
n_iter = 3000;
N_projections = 200; 
GPU = 1; 
cameras = zeros(N_projections,3);
cameras(:,2)=(pi/180)*(0:360/N_projections:360-360/N_projections);

% Simulate projection  
fprintf('Simulating projection ...\n');
projection = et_project(0, cameras, attenuation_phantom, 0, GPU,1);

% Reconstruct volume 
error = zeros(1,n_iter); 
log_likelihood = zeros(1,n_iter); 
fprintf('MLEM reconstruction ...\n');
mask = ones(Nx,Ny,Nz);
B = et_backproject(projection, cameras, 0, 0, GPU); 
attenuation = 0.01*ones(Nx,Ny,Nz); 
for i =1:n_iter
    fprintf('iter %d \n',i);
    update = (et_backproject(et_project(0, cameras, attenuation, 0, GPU,1), cameras, 0, 0, GPU) + 0.0001) ./ (B + 0.0001);
    attenuation = mask.*(attenuation.*update);
    imagesc(squeeze(attenuation(:,Ny/2,:))); axis equal tight off; pause(0.1); 
    error(i) = norm(attenuation(:)-attenuation_phantom(:)); 
end
fprintf('tt_demo_parallel_large done.\n');

