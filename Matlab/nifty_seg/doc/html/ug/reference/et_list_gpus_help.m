%% et_list_gpus
% List installed CUDA compatible GPUs.
%
%% Description
% |gpus = et_list_gpus()|
%
% |gpus| is a matrix of size |[n,5]| where each line represents one detected compatible GPU.
% For each GPU the function returns the following information:
% DEVICE ID, COMPUTING CAPABILITY (Gflops), NUMBER OF MULTIPROCESSORS, CLOCK (Hz), GLOBAL MEMORY (Mbytes)
% If no GPUs are found, the function returns 0.
%
%
%% Example
gpus = et_list_gpus();

%% See also
% <et_set_gpu_help.html |et_set_gpu|>,
% <et_rotate_help.html |et_rotate|>, <et_convolve_help.html |et_convolve|>, <et_affine_help.html |et_affine|>,
% <et_project_help.html |et_project|>, <et_backproject_help.html |et_backproject|>,
%
%
%% Source
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC-UCL.
% Gower Street, London, UK

