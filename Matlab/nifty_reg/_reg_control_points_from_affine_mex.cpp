/*
 *  _reg_control_points_from_affine_mex.cpp
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

   /* Check for correct number of arguments. */
   if (!(nrhs==3 || nrhs==4)){
       mexErrMsgTxt("3 or 4 inputs required: Affine transformation [4,4], image size [1,3], control points grid pacing [1,3], [use_gpu]");
   }
   if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
       mexErrMsgTxt("Affine information must be in double floating point precision.");
   if (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
       mexErrMsgTxt("Control Points must be in double floating point precision.");
   if (!(mxGetClassID(prhs[1]) == mxINT32_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Image Size must be INT32 or double");

   const mxClassID cid_affine      = mxGetClassID(prhs[0]);
   const int       dim_affine      = mxGetNumberOfDimensions(prhs[0]); 

   const mxClassID cid_image_size  = mxGetClassID(prhs[1]);
   const int       dim_image_size  = mxGetNumberOfDimensions(prhs[1]); 

   const mxClassID cid_cp_spacing  = mxGetClassID(prhs[2]);
   const int       dim_cp_spacing  = mxGetNumberOfDimensions(prhs[2]); 

   //extract pointers to input matrices
   double* affineTransformation = (double *) (mxGetData(prhs[0]));
   double* gridSpacing          = (double *) (mxGetData(prhs[2]));

   int imageSize_i[3]; 
   if (mxGetClassID(prhs[1]) == mxINT32_CLASS)
       {
       imageSize_i[0] = ( (int *) mxGetData(prhs[1]) )[0]; imageSize_i[1] = ( (int *) mxGetData(prhs[1]) )[1]; imageSize_i[2] = ( (int *) mxGetData(prhs[1]) )[2];
       }
   else
       {
       double* imageSize            = (double *) (mxGetData(prhs[1]));       
       imageSize_i[0] = imageSize[0]; imageSize_i[1] = imageSize[1]; imageSize_i[2] = imageSize[2];
       }

   //convert data types
   float affineTransformation_f[16]; for (int i=0; i<16; i++) affineTransformation_f[i] = affineTransformation[i];
   float gridSpacing_f[3]; gridSpacing_f[0] = gridSpacing[0]; gridSpacing_f[1] = gridSpacing[1]; gridSpacing_f[2] = gridSpacing[2];

   //calculate size of control points grid, in order to allocate output
   int cp_size[3];
   reg_array_compute_control_point_size(cp_size, imageSize_i, gridSpacing_f);

   //allocate output matrix for control points   
   mwSize mw_cp_size[4];
   mw_cp_size[0] = (mwSize)cp_size[0];
   mw_cp_size[1] = (mwSize)cp_size[1];
   mw_cp_size[2] = (mwSize)cp_size[2];
   mw_cp_size[3] = 3;

   plhs[0] =  mxCreateNumericArray(4, mw_cp_size, mxSINGLE_CLASS, mxREAL);
   float *cp_ptr = (float *)(mxGetData(plhs[0]));

   /* Compute control points from affine */
   status = reg_array_bspline_initialiseControlPointGridWithAffine(cp_ptr, affineTransformation_f, imageSize_i, gridSpacing_f);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while initialising control points.");
   return;
}


