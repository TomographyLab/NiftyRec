/*
 *  _seg_initialise.cpp
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
   if (!(nrhs==3)){
      mexErrMsgTxt("3 input required: image filename, mask filename, number_of_classes");
   }
   /* Call NiftySeg function */
//   int size_x = mxGetDimensions(prhs[0])[0]; 
//   int size_y = mxGetDimensions(prhs[0])[1]; 
//   int size_z = mxGetDimensions(prhs[0])[2]; 
//   int n_dimensions = mxGetDimensions(prhs[0])[0]; 
//   PrecisionTYPE* image_data_ptr = (PrecisionTYPE*) mxGetData(prhs[0]);

   char* image_filename = (char*) mxArrayToString(prhs[0]);
   char* mask_filename = (char*) mxArrayToString(prhs[1]);
   int n_classes = (int) (mxGetScalar(prhs[2]));

//   status = seg_array_initialise(n_classes, size_x, size_y, size_z, n_dimensions, image_data_ptr);
   status = seg_array_initialise_fromfile(n_classes, image_filename, mask_filename);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while initialising NiftySeg.");
   return;
}


