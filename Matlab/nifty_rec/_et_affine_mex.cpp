/*
 *  _et_affine_mex.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>
#include <mex.h>

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (!(nrhs==2 || nrhs==3 || nrhs==4)){
      mexErrMsgTxt("2, 3 or 4 inputs required: InputImage, AffineMatrix, EnableGPU, Background");
   }

   const mxClassID cid_image  = mxGetClassID(prhs[0]);
   const int       dim_image  = mxGetNumberOfDimensions(prhs[0]); 

   const mxClassID cid_affine   = mxGetClassID(prhs[1]);
   const int       dim_affine   = mxGetNumberOfDimensions(prhs[1]); 

   int image_size[3];  //size of input image matrix.
   int affine_size[2];       //size of affine matrix.
   float background;      // Background value, defaults to 0. (It is used to fill empty areas when images are rotated)
   int enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.

   int status = 1;        // Return 0 if everython ok, 1 if errors.

   /* The inputs must be noncomplex floating point matrices */
   for (int n=0; n<2; n++)
      {
      if ( mxGetClassID(prhs[n]) != mxDOUBLE_CLASS )
         mexErrMsgTxt("Input must be noncomplex double floating points.");
      }
   /* Check if size of affine matrix is correct */
   affine_size[0] = mxGetDimensions(prhs[1])[0];
   affine_size[1] = mxGetDimensions(prhs[1])[1];
   if (!(affine_size[0] == 4 && affine_size[1] == 4))
      mexErrMsgTxt("Affine matrix must be of size [4 x 4]");        
   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: TransformedImage");
   }          
   /* Check consistency of input */
   switch(dim_image)
      {
      /* Check consistency of input if 2D */
      case 2:      
           image_size[0] = mxGetDimensions(prhs[0])[0];
           image_size[1] = mxGetDimensions(prhs[0])[1];
           image_size[2] = 1;
           if (image_size[0]<2 || image_size[1]<2)
               mexErrMsgTxt("Size of image matrix must be greater then 2");
           break;
      /* Check consistency of input if 3D (and create size of sinogram for return */
      case 3:        
           image_size[0] = mxGetDimensions(prhs[0])[0];
           image_size[1] = mxGetDimensions(prhs[0])[1];
           image_size[2] = mxGetDimensions(prhs[0])[2];
           if (image_size[0]<2 || image_size[1]<2 || image_size[2]<2)
               mexErrMsgTxt("Size of image matrix must be greater then 2");           
           break;        
      default:
           mexErrMsgTxt("image must be either 2D or 3D.");
           break;
      }

   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>=3)
       enable_gpu = (int) (mxGetScalar(prhs[2]));
  
   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=4)
       background = (float) (mxGetScalar(prhs[3]));

   /* Extract pointers to input matrices */
   double *image_ptr = (double *) (mxGetData(prhs[0]));
   double *affine_ptr = (double *) (mxGetData(prhs[1]));

   /* Allocate out matrix */   
   int dim_out;   
   mwSize mw_out_size[3];
   dim_out = dim_image;
   mw_out_size[0] = (mwSize)image_size[0];
   mw_out_size[1] = (mwSize)image_size[1]; 
   mw_out_size[2] = (mwSize)image_size[2]; 

   plhs[0] =  mxCreateNumericArray(dim_out, mw_out_size, mxDOUBLE_CLASS, mxREAL);
   double *out_ptr = (double *)(mxGetData(plhs[0]));

   /* Convert to float */
   float *image_float = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
   float *out_float = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
   float *affine_float = (float*) malloc(affine_size[0]*affine_size[1]*sizeof(float));
   for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
       image_float[i] = image_ptr[i];
   for (int i=0; i<affine_size[0]*affine_size[1]; i++)
       affine_float[i] = affine_ptr[i];

   /* Perform Affine Transform */
   status = et_array_affine(image_float, image_size, out_float, image_size, affine_float, affine_size, background, enable_gpu);

   /* Convert to double */
   for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
       out_ptr[i] = out_float[i];

   /* Shutdown */
   free(image_float);
   free(out_float);
   free(affine_float);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while performing Affine Transform.");
   return;
}


