

% TT_DEMO_MLEM_PARALLEL
%     NiftyRec Demo: iterative reconstruction transmission tomography with parallel rays 
%
%See also
%   TT_DEMO_PROJECT_PARALLEL 
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


n_iter = 1000;
N = 128; 
N_projections = 120; 
GPU = 1; 
which_phantom = 0;

mask = et_spherical_phantom(N,N,N,N/2,1,0,(N+1)/2,(N+1)/2,(N+1)/2);

if which_phantom==0
    %0)uniform activity:
    attenuation_128 = et_load_nifti('attenuation_02_128.nii');
    attenuation_128 = double(attenuation_128.img);
    phantom_attenuation = mask.*et_rotate(attenuation_128,[pi/2,0,0],[64.5,64.5,64.5],1,0)*1e-5; 
else which_phantom==1
    %1)spheres:
    phantom_attenuation = 1e-3*mask .* (et_spherical_phantom(N,N,N, N/8  ,30 ,0  ,N/2   ,N/2,N/4      ) + ...
                                        et_spherical_phantom(N,N,N, N/12 ,20 ,0  ,N/2   ,N/2,1.1*3*N/4) + ...
                                        et_spherical_phantom(N,N,N, N/10 ,40 ,0  ,N/2   ,N/2,1.1*N/2  ) + ... 
                                        et_spherical_phantom(N,N,N, N/10 ,3  ,0 ,3*N/4 ,N/2,1.0*N/2  ) ) ;
    phantom_attenuation = phantom_attenuation.*et_spherical_phantom(N,N,N, N/10 ,0  ,1  ,N/4   ,N/2,1.0*N/2  ); 
end

cameras = zeros(N_projections,3);
cameras(:,2)=(pi/180)*(0:360/N_projections:360-360/N_projections);

projection = et_project(0, cameras, phantom_attenuation, 0, GPU,1);
min_invprojection = min(1./projection(:)); 
max_invprojection = max(1./projection(:)); 
image((-min_invprojection+1./projection(:,:,1))/max_invprojection*64); colormap gray; axis tight off equal; 
for i=1:N_projections
    image((-min_invprojection+1./projection(:,:,i))/max_invprojection*64); colormap gray; axis tight off equal; pause(0.1)
end

B = et_backproject(projection, cameras, 0, 0, GPU); 
attenuation = 0.01*ones(N,N,N); 
for i =1:n_iter
    fprintf('iter %d \n',i);
    update = (et_backproject(et_project(0, cameras, attenuation, 0, GPU,1), cameras, 0, 0, GPU) + 0.0001) ./ (B + 0.0001);
    attenuation = mask.*(attenuation.*update);
    imagesc(squeeze(attenuation(:,N/2,:))); axis equal tight off; pause(0.1); 
end
    
