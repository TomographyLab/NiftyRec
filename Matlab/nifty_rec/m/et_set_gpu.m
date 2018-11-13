function info = et_set_gpu(device_id)
%ET_SET_GPU
%    Set GPU device to be used by et_.. functions.
%
%Description
%    INFO = ET_SET_GPU(DEVICE_ID)
%
%    DEVICE_ID identifies the device. Use ET_LIST_GPUS to obtain a list of 
%    available CUDA compatible GPUs and their DEVICE_ID. 
%
%    The function returns information about the device if selection is successful, 
%    0 otherwise.
%
%Example
%    gpus = et_list_gpus();
%    et_set_gpu(gpus(1,1));
%
%
%See also
%   ET_LIST_GPUS, ET_PROJECT, ET_BACKPROJECT, ET_ROTATE, ET_AFFINE, ET_CONVOLVE
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK

info = et_set_gpu_mex(device_id);

