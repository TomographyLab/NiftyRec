/*
 *  _seg_set_MRF_strength.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */



#include <mex.h>
#include "_seg_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int status;
   double strength;

   if (nrhs!=1)
        mexErrMsgTxt("One scalar (double) input required: MRF_strength");

   strength = (double) mxGetScalar(prhs[0]);

   status = seg_array_set_MRF_strength(strength);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while setting MRF strength for NiftySeg.");
   return;
}

