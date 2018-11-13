function gpus = et_list_gpus()
%ET_LIST_GPUS
%    List installed CUDA compatible GPUs.
%
%Description
%    GPUS = ET_LIST_GPUS()
%
%    GPUS is a matrix of size [n,5] where each line represents one detected compatible GPU.
%    For each GPU the function returns the following information:
%    DEVICE ID, COMPUTING CAPABILITY (Gflops), NUMBER OF MULTIPROCESSORS, CLOCK (Hz), GLOBAL MEMORY (Mbytes)
%    If no GPUs are found, the function returns 0.
%
%
%Example
%    gpus = et_list_gpus();
%
%
%See also
%   ET_SET_GPU
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK

gpus = et_list_gpus_mex();

