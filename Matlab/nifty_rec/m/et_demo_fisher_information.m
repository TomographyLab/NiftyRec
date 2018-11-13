
% ET_DEMO_FISHER_INFORMATION 
%     NiftyRec SPECT demo: compute the Fisher Information matrix on a grid.  
%
%See also
%   ET_DEMO_FISHER_INFORMATION_2D, ET_DEMO_NOISE_INSTANCES, ET_DEMO_MLEM, ET_DEMO_OSEM 
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


%% Test Fisher Information on 3D grid, sphere in uniform background

N              = 128;
N_cameras      = 180;
%arc            = 2*pi;
psfs           = make_3s_psfs(N, 2.4 , 60*2.4, 3.0, 0.01);
xyz_roi        = [30,40,50];
radius_roi     = 15; 
activity_roi   = 100;
activity_bg    = 10;
mask           = et_spherical_phantom(N,N,N,N/2,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
regularisation = 0; % try with a value around 1000
%grid_spacing   = 4;

attenuation = mask.*zeros(N,N,N);


%% 3D Fisher Matrix
cam         = linspace(0,2*pi,N_cameras);
cameras     = [cam',0*cam',0*cam'];

activity    = mask.*et_spherical_phantom(N,N,N,radius_roi,activity_roi,activity_bg,xyz_roi(1),xyz_roi(2),xyz_roi(3)); 
ytrue = et_project(activity,cameras,attenuation,psfs,1);
ytrue(ytrue<=0.001)=0.001;
iytrue = 1./ytrue;

grid_spacing = 8;
grid = zeros(N,N,N);
counter=0;
for k=1:grid_spacing:N
    for i=1:grid_spacing:N
        for j=1:grid_spacing:N
            counter=counter+1;
            grid(k,j,i)=counter;
        end
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1))); 
var = reshape(diag(cov),int32((size(fisher,1))^(1/3)),int32((size(fisher,1))^(1/3)),int32((size(fisher,1))^(1/3))); 
figure;
for i=1:(size(fisher,1))^(1/3)
    subplot(abs(sqrt((size(fisher,1))^(1/3))),abs(sqrt((size(fisher,1))^(1/3))),i); 
    image(10000*var(:,:,i)); axis square off; 
end

et_reset_gpu();

