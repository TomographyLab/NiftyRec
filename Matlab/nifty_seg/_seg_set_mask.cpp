/*
 *  _seg_set_mask.cpp
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

   /* Check for proper number of arguments. */
   if (!(nrhs==1)){
      mexErrMsgTxt("1 input required: mask filename");
   }

   char* filename = (char*) mxArrayToString(prhs[0]);

   status = seg_array_set_mask_fromfile(filename);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while loading NiftySeg segmentation mask.");
   return;
}


