/*
 *  _seg_get_likelihood.cpp
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
   double loglikelihood=0;

   status = seg_array_get_loglikelihood(&loglikelihood);

   plhs[0] = mxCreateDoubleScalar(loglikelihood);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while getting log likelihood.");
   return;
}

