/*
 *  _et_convolve_mex.cpp
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
   if (!(nrhs==2 || nrhs==3)){
      mexErrMsgTxt("2 or 3 inputs required: Image, PointSpreadFunction, EnableGPU");
   }

   const mxClassID cid_image  = mxGetClassID(prhs[0]);
   const int       dim_image  = mxGetNumberOfDimensions(prhs[0]); 

   const mxClassID cid_psf       = mxGetClassID(prhs[1]);
   const int       dim_psf       = mxGetNumberOfDimensions(prhs[1]); 

   int image_size[3];     //size of input image matrix.
   int psf_size[3];       //size of input psf matrix.
   int out_size[3];       //size of output image matrix.
   int enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.

   int N;                 // Size of image: [N,N,m].
   int m;                 // Size of image: [N,N,m].
   int status = 1;        // Return 0 if everython ok, 1 if errors.

   /* The inputs must be noncomplex floating point matrices */
   for (int n=0; n<2; n++)
      {
      if ( mxGetClassID(prhs[n]) != mxDOUBLE_CLASS )
         mexErrMsgTxt("Input must be noncomplex double floating point.");
      }    
   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: out");
   }          
   /* Check consistency of input (and create size of out for return */
   switch(dim_image)
      {
      /* Check consistency of input if 2D (and create size of out for return) */
      case 2:
           if (dim_psf != 2)
               mexErrMsgTxt("Size of Point Spread Function matrix must match the size of image (ddpsf).");        
           image_size[0] = mxGetDimensions(prhs[0])[0];
           image_size[1] = mxGetDimensions(prhs[0])[1];
           image_size[2] = 1;
           psf_size[0] = mxGetDimensions(prhs[1])[0];
           psf_size[1] = mxGetDimensions(prhs[1])[1];
           psf_size[2] = 1;
           N = image_size[0];
           if (N<2)
               mexErrMsgTxt("Size of image matrix must be greater then 2");
           if (image_size[1]!=N)
               mexErrMsgTxt("2D image must be square (3D image can be of size NxNxm)");
           if ((psf_size[0]%2!=1) || (psf_size[1]%2!=1))
                mexErrMsgTxt("Point Spread Function must be of size hxk; h,k odd");
           out_size[0] = N;
           out_size[1] = N;
           out_size[2] = 1;
           break;
      /* Check consistency of input if 3D (and create size of out for return */
      case 3:
           if (dim_psf != 3)
               mexErrMsgTxt("Size of Point Spread Function matrix must match the size of image (ddpsf).");           
           image_size[0] = mxGetDimensions(prhs[0])[0];
           image_size[1] = mxGetDimensions(prhs[0])[1];
           image_size[2] = mxGetDimensions(prhs[0])[2];
           psf_size[0] = mxGetDimensions(prhs[1])[0];
           psf_size[1] = mxGetDimensions(prhs[1])[1];
           psf_size[2] = mxGetDimensions(prhs[1])[2];
           N = image_size[0];
           m = image_size[2];
           if (N<2 || m<2)
               mexErrMsgTxt("Size of image matrix must be greater then 2");           
           if (image_size[1]!=N)
               mexErrMsgTxt("3D image must be of size [N,N,m]");
           if ((psf_size[0]%2!=1) || ((psf_size[1]%2!=1)) || psf_size[2]!=m)
               mexErrMsgTxt("Point Spread Function must be of size hxkxm for image of size NxNxm; h,k odd");
           out_size[0] = N;
           out_size[1] = N;
           out_size[2] = m;
           break;        
      default:
           mexErrMsgTxt("Image must be either 2D or 3D.");
           break;
      }

   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>=3)
       enable_gpu = (int) (mxGetScalar(prhs[2]));

   /* Extract pointers to input matrices */
   double *image_ptr = (double *) (mxGetData(prhs[0]));
   double *psf_ptr = (double *) (mxGetData(prhs[1]));

   /* Allocate out matrix */   
   int dim_out;   
   mwSize mw_out_size[3];
   
   dim_out = dim_image;
   mw_out_size[0] = (mwSize)out_size[0];
   mw_out_size[1] = (mwSize)out_size[1]; 
   mw_out_size[2] = (mwSize)out_size[2]; 

   plhs[0] =  mxCreateNumericArray(dim_out, mw_out_size, mxDOUBLE_CLASS, mxREAL);
   double *out_ptr = (double *)(mxGetData(plhs[0]));

   /* Convert to float */
   float *image_float = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
   float *out_float = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));
   float *psf_float = (float*) malloc(psf_size[0]*psf_size[1]*psf_size[2]*sizeof(float));
   for (int i=0; i<image_size[0]*image_size[1]*image_size[2]; i++)
       image_float[i] = image_ptr[i];
   for (int i=0; i<psf_size[0]*psf_size[1]*psf_size[2]; i++)
       psf_float[i] = psf_ptr[i];

   /* Perform convolution */
   status = et_array_convolve(image_float, image_size, out_float, out_size, psf_float, psf_size, enable_gpu);

   /* Convert to double */
   for (int i=0; i<out_size[0]*out_size[1]*out_size[2]; i++)
       out_ptr[i] = out_float[i];

   /* Shutdown */
   free(image_float);
   free(out_float);
   free(psf_float);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while performing Convolution.");
   return;
}


