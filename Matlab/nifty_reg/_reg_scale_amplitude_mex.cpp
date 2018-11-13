/*
 *  _reg_scale_amplitude_mex.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include <mex.h>
#include "_reg_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (!(nrhs==3 || nrhs==4)){
      mexErrMsgTxt("3 or 4 inputs required: Image [N,M,K], Min Value, Max Value, [use_gpu]");
   }

   int   image_size[3];     //
   float min_value, max_value;
   int   enable_gpu = 0;    // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int   status = 1;        // Return 0 if everython ok, 1 if errors.

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   /* The inputs must be noncomplex floating point matrices */
   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Input image must be noncomplex single or double floating point.");
   if (!(mxGetClassID(prhs[1]) == mxSINGLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Spline nodes must be noncomplex single or double floating point.");
   if (!(mxGetClassID(prhs[2]) == mxSINGLE_CLASS || mxGetClassID(prhs[2]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Spline nodes must be noncomplex single or double floating point.");


   /* Check if EnableGPU is specified */
   if (nrhs>=4)
       enable_gpu =  (int) (mxGetScalar(prhs[3]));

   /* Extract pointers to input matrices */

   float *image_ptr; 
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       image_ptr = (float *) (mxGetData(prhs[0]));
   else
       {
       double *image_ptr_d = (double *) (mxGetData(prhs[0]));
       image_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
           image_ptr[i] = image_ptr_d[i];
       }

   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       min_value = ((float*)mxGetData(prhs[1]))[0];
   else
       min_value = ((double*)mxGetData(prhs[1]))[0];

   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       max_value = ((float*)mxGetData(prhs[2]))[0];
   else
       max_value = ((double*)mxGetData(prhs[2]))[0];

   /* Allocate out matrix */   
   mwSize mw_outimage_size[3];
   mw_outimage_size[0] = (mwSize)image_size[0];
   mw_outimage_size[1] = (mwSize)image_size[1];
   mw_outimage_size[2] = (mwSize)image_size[2];

   plhs[0] =  mxCreateNumericArray(3, mw_outimage_size, mxSINGLE_CLASS, mxREAL);
   float *outimage_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute NMI gradient with respect of position of nodes of the splines */

   memcpy((void*)outimage_ptr, (void*)image_ptr, image_size[0]*image_size[1]*image_size[2]*sizeof(float));
   status = reg_array_scale_amplitude(outimage_ptr, image_size, min_value, max_value, enable_gpu);

   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(image_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while scaling the image.");
   return;
}


