/*
 *  _et_joint_histogram_mex.cpp
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

int _joint_histogram(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /////////////////////////////////////////
   // Parse the matlab inputs
   /////////////////////////////////////////
	
   // Pointer to matrix_A
   float *matrix_A = (float *) (mxGetData(prhs[0]));
   // Pointer to matrix_B
   float *matrix_B = (float *) (mxGetData(prhs[1]));

   int status = 1;
   int matrix_dim = mxGetNumberOfDimensions(prhs[0]); 
   int *matrix_size = (int*) malloc(matrix_dim*sizeof(int));
   for (int i=0; i<matrix_dim; i++)
       matrix_size[i] = mxGetDimensions(prhs[0])[i];
   int histogram_size = (int)((float *)(mxGetData(prhs[2])))[0];
   float min_A = ((float *)(mxGetData(prhs[3])))[0];
   float max_A = ((float *)(mxGetData(prhs[3])))[1];
   float min_B = ((float *)(mxGetData(prhs[4])))[0];
   float max_B = ((float *)(mxGetData(prhs[4])))[1];

   /////////////////////////////////////////
   // Declare and allocate the outputs
   /////////////////////////////////////////
   mwSize HistSize[2];
   HistSize[0] = histogram_size;
   HistSize[1] = histogram_size;
   plhs[0] =  mxCreateNumericArray(2, HistSize, mxINT32_CLASS, mxREAL);
   int *joint_histogram = (int *)(mxGetData(plhs[0]));

   // Compute joint histogram, no GPU
   status = et_array_joint_histogram(matrix_A, matrix_B, joint_histogram, matrix_dim, matrix_size, histogram_size, min_A, max_A, min_B, max_B, 0);
   free(matrix_size);
   
   return status;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=5 ){
      mexErrMsgTxt("5 inputs required: matrix_A, matrix_B, n_bins, [min_A, max_A], [min_B, max_B]");
   }

   const mxClassID cid_matrix_A  = mxGetClassID(prhs[0]);
   const int       dim_matrix_A  = mxGetNumberOfDimensions(prhs[0]); 
   const mxClassID cid_matrix_B  = mxGetClassID(prhs[1]);
   const int       dim_matrix_B  = mxGetNumberOfDimensions(prhs[1]);
   int status = 1;

   /* The inputs must be noncomplex single matrices. */
   if (cid_matrix_A != mxSINGLE_CLASS)
       mexErrMsgTxt("Input must be single.");
   if (cid_matrix_B != mxSINGLE_CLASS)
       mexErrMsgTxt("Input must be single.");

   if (dim_matrix_A != dim_matrix_B)
       mexErrMsgTxt("The two input matrices must have same dimensions.");
   
   for (int i=0; i<dim_matrix_A; i++)
       if (mxGetDimensions(prhs[0])[i] != mxGetDimensions(prhs[1])[i])
           mexErrMsgTxt("The two input matrices must have same size.");

   if (mxGetDimensions(prhs[2])[0] != 1 || mxGetDimensions(prhs[2])[1] != 1)
       mexErrMsgTxt("Third parameter must be a scalar: n_bins");
       
   if (mxGetDimensions(prhs[3])[0] != 1 || mxGetDimensions(prhs[3])[1] != 2)
       mexErrMsgTxt("Fourth parameter must be [min_A, max_A]"); 

   if (mxGetDimensions(prhs[4])[0] != 1 || mxGetDimensions(prhs[4])[1] != 2)
       mexErrMsgTxt("Fifth parameter must be [min_B, max_B]"); 

   if (nlhs != 1){
      mexErrMsgTxt("One output: Joint histogram");
   }

   status = _joint_histogram(nlhs, plhs, nrhs, prhs);
   
   if (status != 0)
   	mexErrMsgTxt("Error while computing joint histogram.");

   return;
}


