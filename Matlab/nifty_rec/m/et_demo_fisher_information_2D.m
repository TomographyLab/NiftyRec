
% ET_DEMO_FISHER_INFORMATION_2D 
%     NiftyRec SPECT demo: compute the Fisher Information matrix on a 2D grid for different grid strides. 
%
%See also
%   ET_DEMO_FISHER_INFORMATION, ET_DEMO_NOISE_INSTANCES, ET_DEMO_MLEM, ET_DEMO_OSEM
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK



%% Test Fisher Information on 2D grid, uniform background vs sphere in
%% uniform background, pi and 2pi camera rotation

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
%activity    = mask.*et_spherical_phantom(N,N,N,radius_roi,activity_roi,activity_bg,xyz_roi(1),xyz_roi(2),xyz_roi(3));
%cam         = linspace(0,arc,N_cameras);
%cameras     = [cam',0*cam',0*cam'];

%grid        = zeros(N,N,N);
%k=xyz_roi(1);
%counter=0;
%for i=1:grid_spacing:N
%    for j=1:grid_spacing:N
%        counter=counter+1;
%        grid(k,j,i)=counter;
%    end
%end

%ytrue = et_project(activity,cameras,attenuation,psfs,1);
%ytrue(ytrue<=0.001)=0.001;
%iytrue = 1./ytrue;
%fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
%cov = inv(fisher+regularisation*eye(size(fisher,1))); 
%var = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
%subplot(3,3,1); image(30000*var);

figure;


%% Without object - 2pi
cam         = linspace(0,2*pi,N_cameras);
cameras     = [cam',0*cam',0*cam'];

activity    = mask.*et_spherical_phantom(N,N,N,radius_roi,activity_bg,activity_bg,xyz_roi(1),xyz_roi(2),xyz_roi(3)); 
ytrue = et_project(activity,cameras,attenuation,psfs,1);
ytrue(ytrue<=0.001)=0.001;
iytrue = 1./ytrue;

% grid spacing 4
grid_spacing = 4;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1))); 
var1 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,1); image(sqrt(4)*25000*var1); axis square off; pause(0.5)

% grid spacing 3
grid_spacing = 3;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var2 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,2); image(sqrt(3)*25000*var2); axis square off; pause(0.5);

% grid spacing 2
grid_spacing = 2;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var3 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,3); image(sqrt(2)*25000*var3); axis square off; pause(0.5);


%% With object - 2pi
cam         = linspace(0,2*pi,N_cameras);
cameras     = [cam',0*cam',0*cam'];

activity    = mask.*et_spherical_phantom(N,N,N,radius_roi,activity_roi,activity_bg,xyz_roi(1),xyz_roi(2),xyz_roi(3)); 
ytrue = et_project(activity,cameras,attenuation,psfs,1);
ytrue(ytrue<=0.001)=0.001;
iytrue = 1./ytrue;

% grid spacing 4
grid_spacing = 4;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var4 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,4); image(sqrt(4)*25000*var4); axis square off; pause(0.5)

% grid spacing 3
grid_spacing = 3;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var5 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,5); image(sqrt(3)*25000*var5); axis square off; pause(0.5)

% grid spacing 2
grid_spacing = 2;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var6 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,6); image(sqrt(2)*25000*var6); axis square off; pause(0.5)


%% With object - pi
cam         = linspace(0,pi,N_cameras);
cameras     = [cam',0*cam',0*cam'];

activity    = mask.*et_spherical_phantom(N,N,N,radius_roi,activity_roi,activity_bg,xyz_roi(1),xyz_roi(2),xyz_roi(3)); 
ytrue = et_project(activity,cameras,attenuation,psfs,1);
ytrue(ytrue<=0.001)=0.001;
iytrue = 1./ytrue;

% grid spacing 4
grid_spacing = 4;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var7 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,7); image(sqrt(4)*25000*var7); axis square off; pause(0.5)

% grid spacing 3
grid_spacing = 3;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var8 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,8); image(sqrt(3)*25000*var8); axis square off; pause(0.5)

% grid spacing 2
grid_spacing = 2;
grid = zeros(N,N,N);
k=xyz_roi(1);
counter=0;
for i=1:grid_spacing:N
    for j=1:grid_spacing:N
        counter=counter+1;
        grid(k,j,i)=counter;
    end
end
fisher = et_fisher_grid_invprojection_mex(100000*iytrue,cameras,grid,attenuation,psfs,1,0.0001);
cov = inv(fisher+regularisation*eye(size(fisher,1)));  
var9 = reshape(diag(cov),sqrt(size(fisher,1)),sqrt(size(fisher,1))); 
subplot(3,3,9); image(sqrt(2)*25000*var9); axis square off; pause(0.5)

et_reset_gpu();


