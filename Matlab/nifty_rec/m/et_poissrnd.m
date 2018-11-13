function k = et_poissrnd(lambda)
% ET_POISSRND
%     NiftyRec: equivalent to Matworks function poissrnd.m, but does not depend on the 
%     Statistical Toolbox. Wraps Fessler's IRT poisson.m function
%
%See also
%     poissrnd
%
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK

k = poisson(lambda);

