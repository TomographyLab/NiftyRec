/*
 *  _seg_set_biasfield_parameters.cpp
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
   int order;
   double threshold;

   if (nrhs!=2)
        mexErrMsgTxt("Two scalar (double) input parameters required: polynomial_order, threshold");

   order = (int) mxGetScalar(prhs[0]);
   threshold = (double) mxGetScalar(prhs[1]);

   status = seg_array_set_biasfield_parameters(order, threshold);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while setting MRF strength for NiftySeg.");
   return;
}

