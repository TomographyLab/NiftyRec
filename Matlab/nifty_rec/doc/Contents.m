% NiftyRec - Emission Tomography Toolbox
% Version 13-Mar-2012
%
%
% Projection and Backprojection.
% et_project               - Projection to detector space
% et_backproject           - Backprojection from detector space
%
% Reconstruction
% et_mapem_step
% et_osmapem_step
% 
% GPU Acceleration.
% et_list_gpus                 - List available CUDA GPU devices and their compute capability 
% et_set_gpu                   - Select CUDA compatible GPU device
% et_rotate                    - Rotate 2D or 3D image
% et_affine                    - Affine Transform of 2D or 3D image
% et_convolve                  - Multiple 1D or 2D convolutions
% et_project                   - Projection to detector space
% et_backproject               - Backprojection from detector space
% et_mapem_step                - Maximum A Posteriori Expectation Maximisation step
% et_osmapem_step              - Ordered Subsets Maximum A Posteriori Expectation Maximisation step
%
% Gaussian Mixture.
% et_em_gaussian_mixture       - Gaussian Mixture segmentation
%
% Joint Entropy.
% et_je_gradient               - Gradient of Joint Entropy similarity measure
%
% Fisher Information
% et_fisher_grid               - Estimate the Fisher Information matrix from the object
% et_fisher_grid_invprojection - Estimate the Fisher Information matrix from the inverse of the projection (faster)
%
% Miscellaneous
% et_rotate                    - Rotate 2D or 3D image
% et_affine                    - Affine Transform of 2D or 3D image
% et_convolve                  - Multiple 1D or 2D convolutions
% et_resize                    - Resize 2D or 3D image
% et_zeropad                   - Zeropad 2D or 3D image
% et_apply_lesions             - Apply lesions to activity image
% et_spherical_phantom         - Create spherical activity phantom 
%
% Demo
% et_mlem_demo                 - MLEM (Maximum Likelihood Expectation Maximisation) reconstruction demo
% et_osem_demo                 - OSEM (Ordered Subsets Expetation Maximisation) reconstruction demo
% et_mapem_scatter_rbsc_demo   - MLEM (Maximum Likelihood Expectation Maximisation) reconstruction demo with scatter correction
%
%
%
%
% Stefano Pedemonte
% Copyright 2009-2012 CMIC, UCL
% Gower Street, London, UK


