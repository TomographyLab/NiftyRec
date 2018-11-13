/*
 *  _seg_erode.cpp
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
   /* Check for proper number of arguments. */
   if (!(nrhs==2)){
      mexErrMsgTxt("2 inputs required: Image, dilation iterations");
   }

   const int       dim_image  = mxGetNumberOfDimensions(prhs[0]); 

   int image_size[3];     // size of input image matrix.
   int status = 1;        // Return 0 if everython ok, 1 if errors.
   int iterations = 0;    // Dilate itearations

   /* Check dimensions of the input image */
   if (dim_image != 3) mexErrMsgTxt("Erosion of 3D images only.");

   /* The inputs must be noncomplex floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'InputImage' must be noncomplex single or double.");

   /* Check if size of affine matrix is correct */
   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: ErodedImage");
   }          

   /* Extract pointers to input matrices */
   float *image_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       image_ptr = (float *) (mxGetData(prhs[0]));
   else
   {
       double *image_ptr_d = (double *) (mxGetData(prhs[0]));
       image_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2];i++)
           image_ptr[i] = image_ptr_d[i];
   }
   iterations = (int) (mxGetScalar(prhs[1])); 

   /* Allocate out matrix */   
   int dim_out;   
   mwSize mw_out_size[3];
   dim_out = dim_image;
   mw_out_size[0] = (mwSize)image_size[0];
   mw_out_size[1] = (mwSize)image_size[1]; 
   mw_out_size[2] = (mwSize)image_size[2]; 
   plhs[0] =  mxCreateNumericArray(dim_out, mw_out_size, mxSINGLE_CLASS, mxREAL);
   float *out_ptr = (float *)(mxGetData(plhs[0]));

   /* Erode */
   status = seg_array_erode(image_ptr, image_size, out_ptr, iterations);

   /* Shutdown */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(image_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while performing Erosion.");
   return;
}


