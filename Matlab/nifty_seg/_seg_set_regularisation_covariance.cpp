/*
 *  _seg_set_regularisation_covariance.cpp
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
   double regularisation;

   if (nrhs!=1)
        mexErrMsgTxt("One scalar (double) input required: regularisation parameter");

   regularisation = (double) mxGetScalar(prhs[0]);

   status = seg_array_set_regularisation_covariance(regularisation);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("NiftySeg: Error while setting the regularisation parameter for the inversion of the covariance matrix.");
   return;
}

