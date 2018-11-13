% NiftySeg - Segmentation Toolbox
% Version 13-Mar-2012
%
%
% Gaussian Mixture
% seg_initialise                     - Initialise segmenter (see seg_terminate)
% seg_set_biasfield_parameters       - Set bias field parameter
% seg_set_input_image                - Set input image to be segmented
% seg_set_mask                       - Set binary mask for the region to be segmented
% seg_set_mean                       - Set mu parameters of the Gaussian Mixture
% seg_set_MRF_strength               - Set Markov Random Field parameter (spatial dependence)
% seg_set_priors                     - Set the multinomial priors for the tissue labels
% seg_set_regularisation_covariance  - Regularisation parameter for the estimation of the classes covariance matrix 
% seg_set_segmentation               - Set the partial membership manually
% seg_set_variance                   - Set parameter sigma of the Gaussian Mixture
% seg_get_biasfield                  - Get bias field corrected image 
% seg_get_biasfield_parameter        - Get bias field correction parameters 
% seg_get_image                      - Get input image
% seg_get_loglikelihood              - Get log likelihood
% seg_get_mask                       - Get binary mask of the region to be segmented 
% seg_get_mean_variance              - Get the parameters of the Gaussian Mixture
% seg_get_MRF_strength               - Get Markov Random Field parameter (spatial dependence)
% seg_get_priors                     - Get multinomial priors (statistical atlas)
% seg_get_regularisation_covariance  - Get regularisation parameter for the estimation of the classes covariance matrix 
% seg_get_segmentation               - Get the partial membership 
% seg_save_biasfield                 - Save to file the bias field corrected image 
% seg_save_mask                      - Save to file the binary mask of the region to be segmented 
% seg_save_priors                    - Save the probabilistic atlas to file
% seg_save_segmentation              - Save the probabilistic segmentation to file
% seg_step                           - Generalised Expectation Maximisation step (includes all steps)
% seg_step_BiasField                 - Update bias field parameters only
% seg_step_Expectation               - Update expectation only
% seg_step_Maximisation              - Update Maximisation only 
% seg_step_Gaussian                  - Update Gaussian Mixture parameters only
% seg_step_MRF                       - Update Markov Random Field only
% seg_step_PriorWeight               - Update prior weight only
% seg_terminate                      - Terminate segmenter (see seg_initialise)
%
% Utilities
% seg_erode                          - Erode segmentation
%
% Demo
% seg_demo_manual_priors             - Generalised Expectation Maximisation segmentation with manual priors 
% seg_demo_auto_priors               - Generalised Expectation Maximisation segmentation with automatic priors 
%
%
% Stefano Pedemonte 
% Copyright 2009-2012 CMIC, UCL
% Gower Street, London, UK


