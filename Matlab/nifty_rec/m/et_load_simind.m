function [sinogram] = et_load_simind(filename,size_proj_X,size_proj_Y,N_projections)
% ET_LOAD_SIMIND
%     NiftyRec: Load Simind simulation output files. 
%
%See also
%     et_load_nifti, et_load_ecat
%
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


fid = fopen(filename);
sinogram = reshape(fread(fid,'float32'),size_proj_X,size_proj_Y,N_projections);
fclose(fid);



