
niftyseg_path = '/Users/pidi/Desktop/niftyinstall/niftyseg/matlab/';
priors_path   = '/Users/pidi/Desktop/niftyinstall/priors/';

addpath(niftyseg_path);

N_iter = 20;
input_image = strcat(priors_path,'T1.nii');
mask_image  = strcat(priors_path,'T1.nii');
out_image   = './niftyseg_segmentation.nii';
priors      = {  strcat(priors_path,'WM_prior.nii'), ...
                 strcat(priors_path,'ECSF_prior.nii'), ...
                 strcat(priors_path,'DGM_prior.nii'), ...
                 strcat(priors_path,'CGM_prior.nii'), ... 
                 strcat(priors_path,'ICSF_prior.nii')};

seg_initialise(input_image, mask_image, 5);
seg_set_MRF_strength(0.5);
seg_set_biasfield_parameters(4,0);
seg_set_priors(priors{1},priors{2},priors{3},priors{4},priors{5});

figure;

for iter = 1:N_iter
    fprintf('Step: %d\n',iter);
    seg_step()
    seg = seg_get_segmentation();
    subplot(2,3,1); imagesc(seg(:,:,50,1)); colormap gray;
    subplot(2,3,2); imagesc(seg(:,:,50,2)); colormap gray;
    subplot(2,3,3); imagesc(seg(:,:,50,3)); colormap gray;
    subplot(2,3,4); imagesc(seg(:,:,50,4)); colormap gray;
    subplot(2,3,5); imagesc(seg(:,:,50,5)); colormap gray;
    pause(0.1);
end

seg_save_segmentation(out_image);



