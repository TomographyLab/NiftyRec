/*
 *  _reg_image_gradient_mex.cpp
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
      mexErrMsgTxt("3 or 4 inputs required: Image [N,M,K], Spline Nodes [Y,W,Z,4], Nodes Spacing [y,w,z], [use_gpu]");
   }

   int   image_size[3];     //
   float spacing[3];        //
   int   nodes_size[3];
   int   enable_gpu = 1;    // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int   status = 1;        // Return 0 if everython ok, 1 if errors.

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   nodes_size[0] = mxGetDimensions(prhs[1])[0];
   nodes_size[1] = mxGetDimensions(prhs[1])[1];
   nodes_size[2] = mxGetDimensions(prhs[1])[2];

   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       {
       spacing[0] = ((float*) mxGetData(prhs[2]))[0];
       spacing[1] = ((float*) mxGetData(prhs[2]))[1];
       spacing[2] = ((float*) mxGetData(prhs[2]))[2];
       }
   else if (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS)
       {
       spacing[0] = ((double*) mxGetData(prhs[2]))[0];
       spacing[1] = ((double*) mxGetData(prhs[2]))[1];
       spacing[2] = ((double*) mxGetData(prhs[2]))[2];
       }

   /* The inputs must be noncomplex floating point matrices */
   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Input image must be noncomplex single or double floating point.");
   if (!(mxGetClassID(prhs[1]) == mxSINGLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
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
   float *nodes_ptr;    
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       nodes_ptr = (float *) (mxGetData(prhs[1]));
   else
       {
       double *nodes_ptr_d = (double *) (mxGetData(prhs[1]));
       nodes_ptr = (float*) malloc(nodes_size[0]*nodes_size[1]*nodes_size[2]*3*sizeof(float));
       for (int i=0; i<nodes_size[0]*nodes_size[1]*nodes_size[2]*3; i++)
           nodes_ptr[i] = nodes_ptr_d[i];
       }

   /* Allocate out matrix */   
   mwSize mw_outimage_size[4];
   mw_outimage_size[0] = (mwSize)image_size[0];
   mw_outimage_size[1] = (mwSize)image_size[1];
   mw_outimage_size[2] = (mwSize)image_size[2];
   mw_outimage_size[3] = 3;

   plhs[0] =  mxCreateNumericArray(4, mw_outimage_size, mxSINGLE_CLASS, mxREAL);
   float *outimage_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute NMI gradient with respect of position of nodes of the splines */

   status = reg_array_image_gradient(outimage_ptr, image_ptr, nodes_ptr, image_size, nodes_size, spacing, enable_gpu);

   if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) free(image_ptr);
   if (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) free(nodes_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing the gradient of the image.");
   return;
}


