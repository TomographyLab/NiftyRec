/*
 *  _seg_get_biasfield_parameters.cpp
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
   double strength=0;

   if (nrhs!=0)
        mexErrMsgTxt("No inputs expected.");

   int enabled; 
   unsigned int order; 
   double threshold;

   status = seg_array_get_biasfield_parameters(&enabled, &order, &threshold);

   plhs[0] = mxCreateDoubleScalar(enabled);
   plhs[1] = mxCreateDoubleScalar(order);
   plhs[2] = mxCreateDoubleScalar(threshold);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while getting bias field correction parameters.");
   return;
}

