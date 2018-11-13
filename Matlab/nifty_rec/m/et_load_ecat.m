function [imaVOL, scaninfo, hd] = et_load_ecat(filename)
% ET_LOAD_ECAT
%     NiftyRec: Load Siemens ECAT files. Wraps loadecat.m, 
%     based on ecatfile.m and readecatvol.m by Flemming
%     Hermansen. 
%
%See also
%     et_load_nifti, et_load_simind
%
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK

[imaVOL,scaninfo,hd]= loadecat(filename);
