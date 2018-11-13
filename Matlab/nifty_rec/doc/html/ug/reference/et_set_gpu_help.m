%% et_set_gpu
% Set GPU device to be used by et_.. functions.
%
%% Description
% |info = et_set_gpu(device_id)|
%
% |device_id| identifies the device. Use |et_list_gpus| to obtain a list of 
% available CUDA compatible GPUs and their |device_id|. 
%
% The function returns information about the device if selection is successful, 
% |0| otherwise.
%
%% Example
gpus = et_list_gpus();
et_set_gpu(gpus(1,1));

%% See also
% <et_list_gpus_help.html |et_list_gpus|>,
% <et_rotate_help.html |et_rotate|>, <et_convolve_help.html |et_convolve|>, <et_affine_help.html |et_affine|>,
% <et_project_help.html |et_project|>, <et_backproject_help.html |et_backproject|>,
%
%
%% Source 
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC-UCL.
% Gower Street, London, UK


