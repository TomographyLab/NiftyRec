/*
 *  _seg_get_regularisation_covariance.cpp
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
   double regularisation=0;

   if (nrhs!=0)
        mexErrMsgTxt("No inputs expected.");

   status = seg_array_get_regularisation_covariance(&regularisation);

   plhs[0] = mxCreateDoubleScalar(regularisation);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt(" NiftySeg: Error while getting the regularisation parameter for the inversion of the covariance matrix.");
   return;
}

