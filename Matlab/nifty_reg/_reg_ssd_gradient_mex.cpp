/*
 *  _reg_ssd_gradient_mex.cpp
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
   int status = 1;
   int use_gpu = 1;
   int image_size[3];

   /* Check for correct number of arguments. */
   if (!(nrhs==4 || nrhs==5)){
      mexErrMsgTxt("4 or 5 inputs required: Target Image, Source Image, Spatial Gradient of the deformed Source Image, Smoothing radius, [use_gpu]");
   }

   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Target Image must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[1]) == mxSINGLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Source Image must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[2]) == mxSINGLE_CLASS || mxGetClassID(prhs[2]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Spatial Gradient of the deformed Source Image must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[3]) == mxINT32_CLASS || mxGetClassID(prhs[2]) == mxSINGLE_CLASS || mxGetClassID(prhs[2]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Smoothing Radius must be real INT32.");

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];


   /* Check if EnableGPU is specified */
   if (nrhs>=5)
       use_gpu =  (int) (mxGetScalar(prhs[4]));

   /* Extract pointers to input matrices */
   float *target_ptr; 
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       target_ptr = (float *) (mxGetData(prhs[0]));
   else
       {
       double *target_ptr_d = (double *) (mxGetData(prhs[0]));
       target_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
           target_ptr[i] = target_ptr_d[i];
       }

   float *source_ptr; 
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       source_ptr = (float *) (mxGetData(prhs[1]));
   else
       {
       double *source_ptr_d = (double *) (mxGetData(prhs[1]));
       source_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
           source_ptr[i] = source_ptr_d[i];
       }

   float *gradient_ptr; 
   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       gradient_ptr = (float *) (mxGetData(prhs[2]));
   else
       {
       double *gradient_ptr_d = (double *) (mxGetData(prhs[2]));
       gradient_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*3*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2]*3; i++)
           gradient_ptr[i] = gradient_ptr_d[i];
       }

   int smoothing_radius[3];
   if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
       {
       smoothing_radius[0] = ((float*) mxGetData(prhs[3]))[0];
       smoothing_radius[1] = ((float*) mxGetData(prhs[3]))[1];
       smoothing_radius[2] = ((float*) mxGetData(prhs[3]))[2];
       }
   else if (mxGetClassID(prhs[3]) == mxDOUBLE_CLASS)
       {
       smoothing_radius[0] = ((double*) mxGetData(prhs[3]))[0];
       smoothing_radius[1] = ((double*) mxGetData(prhs[3]))[1];
       smoothing_radius[2] = ((double*) mxGetData(prhs[3]))[2];
       }
   else if (mxGetClassID(prhs[3]) == mxINT32_CLASS)
       {
       smoothing_radius[0] = ((int*) mxGetData(prhs[3]))[0];
       smoothing_radius[1] = ((int*) mxGetData(prhs[3]))[1];
       smoothing_radius[2] = ((int*) mxGetData(prhs[3]))[2];
       }


   //allocate output matrix for SSD control point gradient
   mwSize mw_size[4];
   mw_size[0] = (mwSize)image_size[0];
   mw_size[1] = (mwSize)image_size[1];
   mw_size[2] = (mwSize)image_size[2];
   mw_size[3] = 3;

   plhs[0] =  mxCreateNumericArray(4, mw_size, mxSINGLE_CLASS, mxREAL);
   float *ssd_gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute the gradient */
   status = reg_array_ssd_gradient(ssd_gradient_ptr, target_ptr, source_ptr, gradient_ptr, image_size, smoothing_radius, use_gpu);

   /* Dealloc */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
       free(target_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
       free(source_ptr);
   if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
       free(gradient_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing SSD gradient wrt control points position.");
   return;
}


