function status = et_isinstalled()
%ET_ISINSTALLED
%    Checks if NiftyRec libraries are installed. 
%
%Description
%    STATUS = ET_ISINSTALLED()
%
%    STATUS is 1 if NiftyRec libraries are installed and functioning, 0 otherwise. 
%
%See also
%   ET_LIST_GPUS
%
% 
%Stefano Pedemonte
%Copyright 2009-2010 CMIC-UCL.
%Gower Street, London, UK


status = et_isinstalled_mex();


