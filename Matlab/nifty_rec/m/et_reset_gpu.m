function status = et_reset_gpu()
%ET_RESET_GPU
%    NiftyRec: reset Graphics Processing Unit device. Call this function before 'clear all' 
%    in case 'clear all' should cause segmentation fault. 
%
%Description
%    STATUS = ET_RESET_GPU()
%
%See also
%   ET_SET_GPU, ET_LIST_GPUS,
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK

status = et_reset_gpu_mex();




