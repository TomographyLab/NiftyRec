/*
 *  _et_rotate_mex.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <mex.h>
#include "_et_array_interface.h"

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
   if (!(nrhs==3 || nrhs==4 || nrhs==5)){
      mexErrMsgTxt("3, 4 or 5 inputs required: InputImage, Rotation, Center, EnableGPU, Background");
   }

   const mxClassID cid_image  = mxGetClassID(prhs[0]);
   const int       dim_image  = mxGetNumberOfDimensions(prhs[0]); 

   const mxClassID cid_rotation   = mxGetClassID(prhs[1]);
   const int       dim_rotation   = mxGetNumberOfDimensions(prhs[1]); 

   const mxClassID cid_center    = mxGetClassID(prhs[2]);
   const int       dim_center    = mxGetNumberOfDimensions(prhs[2]); 

   int image_size[3];     // Size of input image matrix.
   int rotation_size[2];  // Size of rotation vector.
   int center_size[2];    // Size of centers vector.
   float background;      // Background value, defaults to 0. (It is used to fill empty areas when images are rotated).
   int enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.

   int status = 1;        // Return 0 if everython ok, 1 if errors.

   /* The inputs must be noncomplex floating point matrices */
   for (int n=0; n<2; n++)
      {
      if ( mxGetClassID(prhs[n]) != mxDOUBLE_CLASS )
         mexErrMsgTxt("Input must be noncomplex double floating points.");
      }
   /* Check if size of rotation and center is correct */
   rotation_size[0] = mxGetDimensions(prhs[1])[0];
   rotation_size[1] = mxGetDimensions(prhs[1])[1];
   center_size[0] = mxGetDimensions(prhs[2])[0];
   center_size[1] = mxGetDimensions(prhs[2])[1];
   if (rotation_size[1] != dim_image)
      mexErrMsgTxt("Rotation vector must be of length 2 for 2D image and length 3 for 3D image");    
   if (center_size[1] != dim_image)
      mexErrMsgTxt("Center vector must be of length 2 for 2D image and length 3 for 3D image");        
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
   if (nrhs>=4)
       enable_gpu = (int) (mxGetScalar(prhs[3]));
  
   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=5)
       background = (float) (mxGetScalar(prhs[4]));

   /* Extract pointers to input matrices */
   double *image_ptr = (double *) (mxGetData(prhs[0]));
   double *rotation_ptr = (double *) (mxGetData(prhs[1]));
   double *center_ptr = (double *) (mxGetData(prhs[2]));

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
   float *rotation_float = (float*) malloc(rotation_size[0]*rotation_size[1]*sizeof(float));
   float *center_float = (float*) malloc(center_size[0]*center_size[1]*sizeof(float));
   for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
       image_float[i] = image_ptr[i];
   for (int i=0; i<rotation_size[0]*rotation_size[1]; i++)
       rotation_float[i] = -rotation_ptr[i]; //negative
   for (int i=0; i<center_size[0]*center_size[1]; i++)
       center_float[i] = center_ptr[i];

   /* Shift 'center' by 1 (Matlab to C indexes) */
   for (int i=0; i<center_size[0]*center_size[1]; i++)
       center_float[i] = center_float[i]-1;

   /* Perform Affine Transform */
   status = et_array_rotate(image_float, image_size, out_float, rotation_float, center_float, background, ZYX_ROTATION, enable_gpu);
   /* Convert to double */
   for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
       out_ptr[i] = out_float[i];

   /* Shutdown */
   free(image_float);
   free(out_float);
   free(rotation_float);
   free(center_float);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while performing Rotation.");
   return;
}



