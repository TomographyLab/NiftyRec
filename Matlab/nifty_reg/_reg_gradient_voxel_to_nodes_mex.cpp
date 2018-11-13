/*
 *  _reg_gradient_voxel_to_nodes_mex.cpp
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
   int nodes_size[3];

   /* Check for correct number of arguments. */
   if (!(nrhs==3 || nrhs==4)){
      mexErrMsgTxt("3 or 4 inputs required: Gradient wrt voxel intensity [N,M,K], Control points [Y,W,Z,3], Control Points Spacing [y,w,z], [use_gpu]");
   }

   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Gradient must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[1]) == mxSINGLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Spline Nodes must be real single or double precision floating point.");

   if (!(mxGetClassID(prhs[2]) == mxSINGLE_CLASS || mxGetClassID(prhs[2]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Initial Distance Between Nodes must be real single or double precision floating point.");

   image_size[0] = mxGetDimensions(prhs[0])[0];
   image_size[1] = mxGetDimensions(prhs[0])[1];
   image_size[2] = mxGetDimensions(prhs[0])[2];

   nodes_size[0] = mxGetDimensions(prhs[1])[0];
   nodes_size[1] = mxGetDimensions(prhs[1])[1];
   nodes_size[2] = mxGetDimensions(prhs[1])[2];

   /* Check if EnableGPU is specified */
   if (nrhs>=4)
       use_gpu =  (int) (mxGetScalar(prhs[3]));

   /* Extract pointers to input matrices */
   float spacing[3];
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

   float *gradient_ptr; 
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       gradient_ptr = (float *) (mxGetData(prhs[0]));
   else
       {
       double *gradient_ptr_d = (double *) (mxGetData(prhs[0]));
       gradient_ptr = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*3*sizeof(float));
       for (int i=0; i<image_size[0]*image_size[1]*image_size[2]*3; i++)
           gradient_ptr[i] = gradient_ptr_d[i];
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

   int nodes_size_in[3];
   reg_array_compute_control_point_size(nodes_size_in, image_size, spacing);
   if (!(nodes_size_in[0]==nodes_size[0] && nodes_size_in[1]==nodes_size[1] && nodes_size_in[2]==nodes_size[2] ))
       mexErrMsgTxt("Size of nodes is not consistent with nodes spacing. ");

   //allocate output matrix for control points gradient  
   mwSize mw_cp_size[4];
   mw_cp_size[0] = (mwSize)nodes_size_in[0];
   mw_cp_size[1] = (mwSize)nodes_size_in[1];
   mw_cp_size[2] = (mwSize)nodes_size_in[2];
   mw_cp_size[3] = 3;

   plhs[0] =  mxCreateNumericArray(4, mw_cp_size, mxSINGLE_CLASS, mxREAL);
   float *cp_gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute control points from affine */
   status = reg_array_gradient_voxel_to_nodes(cp_gradient_ptr, gradient_ptr, nodes_ptr, image_size, nodes_size, spacing, use_gpu);

   /* Dealloc */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
       free(gradient_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
       free(nodes_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while computing gradient wrt control points position.");
   return;
}


