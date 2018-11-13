function [fwhm slope efficiency] = et_collimator_parallelholes(d,l,h,t,intrinsic,distance,material)
%%% [PSF fwhm slope efficiency] = et_collimator_parallelholes(d,l,h,t,intrinsic,distance,material)

%ET_COLLIMATOR_PARALLELHOLES
%    Creates Point Spread function for parallel holes collimator for SPECT. 
%
%Description
%    Creates Point Spread function for parallel holes collimator for SPECT. 
%
%    [PSF fwhm slope efficiency] = ET_COLLIMATOR_PARALLELHOLES(d,l,h,t,intrinsic,distance)
%
%See also
%    ET_PROJECT, ET_MLEM_DEMO, ET_OSEM_DEMO
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK

if not(exist('material','var'))
    material = 'lead';
end

if material == 'lead'
    mu = 2.17;
end

le = l-2/mu;

fwhm = sqrt((d^2*((distance+le)^2)/le^2) + intrinsic^2); 
slope = sqrt(fwhm^2-intrinsic^2)/distance;
efficiency = d^4/(4*pi*le^2*(d+t)^2);
%PSF = make_3s_psfs(N_y, dx , obj2det, intrinsic, slope);
% ...