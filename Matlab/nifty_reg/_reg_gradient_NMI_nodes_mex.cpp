/*
 *  _reg_gradient_NMI_nodes_mex.cpp
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
   if (!(nrhs==5 || nrhs==6)){
      mexErrMsgTxt("5 or 6 inputs required: Target Image, Source Image, Spline nodes, Initial distance between nodes, Number of histogram bins, [use_gpu]");
   }

   int image_size[3];     //
   float spacing[3];      //
   int nodes_size[4];
   int enable_gpu = 1;    // Flag for GPU Acceleration: 1 to enable, 0 to disable.

   int status = 1;        // Return 0 if everython ok, 1 if errors.

   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Target Image must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[1]) == mxSINGLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Source Image must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[2]) == mxSINGLE_CLASS || mxGetClassID(prhs[2]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Spline Nodes must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[3]) == mxSINGLE_CLASS || mxGetClassID(prhs[3]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Initial Distance Between Nodes must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[4]) == mxINT32_CLASS || mxGetClassID(prhs[4]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Number of histogram bins must be real double precision floating point or INT32.");

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   nodes_size[0] = mxGetDimensions(prhs[2])[0];
   nodes_size[1] = mxGetDimensions(prhs[2])[1];
   nodes_size[2] = mxGetDimensions(prhs[2])[2];
   nodes_size[2] = mxGetDimensions(prhs[2])[3];

   /* Check if EnableGPU is specified */
   if (nrhs>=6)
       enable_gpu =  (int) (mxGetScalar(prhs[5]));

   /* Extract pointers to input matrices */
   if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
       {
       spacing[0] = ((float*) mxGetData(prhs[3]))[0];
       spacing[1] = ((float*) mxGetData(prhs[3]))[1];
       spacing[2] = ((float*) mxGetData(prhs[3]))[2];
       }
   else if (mxGetClassID(prhs[3]) == mxDOUBLE_CLASS)
       {
       spacing[0] = ((double*) mxGetData(prhs[3]))[0];
       spacing[1] = ((double*) mxGetData(prhs[3]))[1];
       spacing[2] = ((double*) mxGetData(prhs[3]))[2];
       }

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

   float *nodes_ptr;    
   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       nodes_ptr = (float *) (mxGetData(prhs[2]));
   else
       {
       double *nodes_ptr_d = (double *) (mxGetData(prhs[2]));
       nodes_ptr = (float*) malloc(nodes_size[0]*nodes_size[1]*nodes_size[2]*sizeof(float));
       for (int i=0; i<nodes_size[0]*nodes_size[1]*nodes_size[2]; i++)
           nodes_ptr[i] = nodes_ptr_d[i];
       }

   int binning = floor(((double*)(mxGetData(prhs[4])))[0]);


fprintf(stderr, "- Nodes size: %d %d %d\n", nodes_size[0], nodes_size[1], nodes_size[2]);
fprintf(stderr, "- Image size: %d %d %d\n", image_size[0], image_size[1], image_size[2]);

   mwSize mw_gradient_size[4];
   mw_gradient_size[0] = 4;
   mw_gradient_size[1] = (mwSize)image_size[0];
   mw_gradient_size[2] = (mwSize)image_size[1];
   mw_gradient_size[3] = (mwSize)image_size[2];

plhs[0] =  mxCreateNumericArray(4, mw_gradient_size, mxSINGLE_CLASS, mxREAL);
float *gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute NMI gradient with respect of position of nodes of the splines */
   status = reg_array_gradient_NMI_nodes(target_ptr ,source_ptr, gradient_ptr, image_size, nodes_ptr, spacing, binning, enable_gpu);

   /* Dealloc */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
       free(target_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
       free(source_ptr);
   if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
       free(nodes_ptr);


   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing NMI gradient with respect of control point positions.");
   return;
}


