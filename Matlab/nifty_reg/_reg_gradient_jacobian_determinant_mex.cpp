/*
 *  _reg_gradient_jacobian_determinant_mex.cpp
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
   int use_gpu = 0;
   int image_size[3];
   int cp_size[3];
   float cp_spacing[3];
   float weight;

   /* Check for correct number of arguments. */
   if (!(nrhs==4 || nrhs==5)){
      mexErrMsgTxt("4 or 5 inputs required: Control points [Y,W,Z,3], Control Points Spacing [1,1,1], Image Size [1,1,1], weight, [use_gpu]");
   }



   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Control Points must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[1]) == mxSINGLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Initial Distance Between Nodes must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[2]) == mxSINGLE_CLASS || mxGetClassID(prhs[2]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Image size must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[3]) == mxSINGLE_CLASS || mxGetClassID(prhs[3]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Weight must be real single or double precision floating point.");

   cp_size[0] = mxGetDimensions(prhs[0])[0];
   cp_size[1] = mxGetDimensions(prhs[0])[1];
   cp_size[2] = mxGetDimensions(prhs[0])[2];

   /* Check if EnableGPU is specified */
   if (nrhs>=5)
       use_gpu =  (int) (mxGetScalar(prhs[4]));

   /* Extract pointers to input matrices */
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       {
       cp_spacing[0] = ((float*) mxGetData(prhs[1]))[0];
       cp_spacing[1] = ((float*) mxGetData(prhs[1]))[1];
       cp_spacing[2] = ((float*) mxGetData(prhs[1]))[2];
       }
   else if (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS)
       {
       cp_spacing[0] = ((double*) mxGetData(prhs[1]))[0];
       cp_spacing[1] = ((double*) mxGetData(prhs[1]))[1];
       cp_spacing[2] = ((double*) mxGetData(prhs[1]))[2];
       }

   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       {
       image_size[0] = ((float*) mxGetData(prhs[2]))[0];
       image_size[1] = ((float*) mxGetData(prhs[2]))[1];
       image_size[2] = ((float*) mxGetData(prhs[2]))[2];
       }
   else if (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS)
       {
       image_size[0] = ((double*) mxGetData(prhs[2]))[0];
       image_size[1] = ((double*) mxGetData(prhs[2]))[1];
       image_size[2] = ((double*) mxGetData(prhs[2]))[2];
       }

   float *control_points_ptr; 
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       control_points_ptr = (float *) (mxGetData(prhs[0]));
   else
       {
       double *control_points_ptr_d = (double *) (mxGetData(prhs[0]));
       control_points_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*3*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2]*3; i++)
           control_points_ptr[i] = control_points_ptr_d[i];
       }

   double weight_d = (double)(mxGetScalar(prhs[3]));
   weight = weight_d;

   int cp_size_in[3];
   reg_array_compute_control_point_size(cp_size_in, image_size, cp_spacing);
   if (!(cp_size_in[0]==cp_size[0] && cp_size_in[1]==cp_size[1] && cp_size_in[2]==cp_size[2] ))
       mexErrMsgTxt("Size of nodes is not consistent with nodes spacing. ");

   //allocate output matrix for control points gradient  
   mwSize mw_cp_size[4];
   mw_cp_size[0] = (mwSize)cp_size_in[0];
   mw_cp_size[1] = (mwSize)cp_size_in[1];
   mw_cp_size[2] = (mwSize)cp_size_in[2];
   mw_cp_size[3] = 3;

   plhs[0] =  mxCreateNumericArray(4, mw_cp_size, mxSINGLE_CLASS, mxREAL);
   float *cp_gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute control points from affine */
   status = reg_array_gradient_jacobian_determinant(cp_gradient_ptr, control_points_ptr, image_size, cp_size, cp_spacing, weight, use_gpu);

   /* Dealloc */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
       free(control_points_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing jacobian determinant gradient wrt control points locations.");
   return;
}


