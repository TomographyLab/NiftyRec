function [nifti_struct] = et_load_nifti(filename)
% ET_LOAD_NIFTI
%     NiftyRec: Load Nifti files. Wraps load_nii.m
%
%See also
%     et_load_ecat, et_load_simind 
%
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK

[nifti_struct]= load_nii(filename);
