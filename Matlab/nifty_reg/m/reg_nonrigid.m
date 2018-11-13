
addpath /home/spedemon/Desktop/niftyinstall/niftyrec/matlab/

N = 192;
cp_distance = [10,10,10]; 

N_iter_per_level = 15;
alpha = 10e-4;
N_levels = 4;
initial_smoothing_sigma = 10;
smoothing_radius = 5;

phantom = et_spherical_phantom(N,N,N,30,100,10,N/2,N/2,N/2); 
control_points = reg_control_points_from_affine_mex(eye(4),[N,N,N],cp_distance); 

control_points_def = control_points; 
%control_points_def(14,11,11,:) = [150,92,92];
control_points_def(:,:,:,1) = control_points(:,:,:,1) + 20*rand(23,23,23);
control_points_def(:,:,:,2) = control_points(:,:,:,2) - 20*rand(23,23,23);
control_points_def(:,:,:,3) = control_points(:,:,:,3) + 20*rand(23,23,23);

target = phantom;
source = reg_resample_spline_mex(target,control_points_def,cp_distance);

control_points_estim = reg_control_points_from_affine_mex(eye(4),[N,N,N],cp_distance); 

smoothing_sigmas = linspace(initial_smoothing_sigma,1,N_levels);
for level = 1:N_levels
    smoothing_sigma = smoothing_sigmas(level);
    fprintf('Level: %2d / %2d   (sigma = %2.2f) \n',level,N_levels,smoothing_sigma);
    target_sm = reg_gaussian_smooth_mex(target,smoothing_sigma);
    source_sm = reg_gaussian_smooth_mex(source,smoothing_sigma);
    result_image = reg_resample_spline_mex(source_sm,control_points_estim,cp_distance);
    alpha_level = alpha / level ;
    for iter = 1:N_iter_per_level
        fprintf('    Iter: %2d / %2d \n',iter,N_iter_per_level);
        grad_deformed = reg_image_gradient_mex(source_sm,control_points_estim,cp_distance,1);
        grad_ssd_voxel = reg_ssd_gradient_mex(target_sm,result_image,grad_deformed,[smoothing_radius,smoothing_radius,smoothing_radius],0);
        grad_ssd_cp = reg_gradient_voxel_to_nodes_mex(grad_ssd_voxel,control_points_estim,cp_distance);
        grad_bending_energy_cp = reg_gradient_bending_energy_mex(control_points_estim,[cp_distance, cp_distance, cp_distance],[N,N,N],5);
        grad_jacobian_det_cp = reg_gradient_jacobian_determinant_mex(control_points_estim,[cp_distance, cp_distance, cp_distance],[N,N,N],0.1);
    
        result_image = reg_resample_spline_mex(source_sm,control_points_estim,cp_distance);
        subplot(2,3,1); imagesc(target_sm(:,:,N/2)); axis tight off; colormap gray; pause(0.05);
        subplot(2,3,2); imagesc(source_sm(:,:,N/2)); axis tight off; colormap gray; pause(0.05);
        subplot(2,3,3); imagesc(result_image(:,:,N/2)); axis tight off; colormap gray; pause(0.05);
        subplot(2,3,4); imagesc(grad_deformed(:,:,N/2,1)); axis tight off; colormap gray; pause(0.05);
        subplot(2,3,5); imagesc(grad_ssd_voxel(:,:,N/2,1)); axis tight off; colormap gray; pause(0.05);
        subplot(2,3,6); imagesc(grad_ssd_cp(:,:,11,1)); axis tight off; colormap gray; pause(0.05);
    
        control_points_estim = control_points_estim - alpha_level * (grad_ssd_cp+grad_bending_energy_cp+grad_jacobian_det_cp); 
    end
end


