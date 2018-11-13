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
/* Matlab extension */
/*###############################################################################################*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Inputs: 
   // 1) input_data: array of input data, single or double floating point. Size: [N].
   // 2) weights: array of weights; weights express the probability that each input value belongs to any of the N_classes. Size: [NxN_classes]. 
   // 3) N_bins: number of bins of the histogram. 
   // 4) min_value: Minimum value for the histogram: all input values that are smaller than min_value are associated to the first histogram bin. 
   // 5) max_value: Maximum value for the histogram: all input values that are larger than max_value are associated to the first histogram bin. 
   // Outputs: 
   // 1) weighted histogram. 

   /* Check correct number of arguments. */
   if (!(nrhs==5)){
      mexErrMsgTxt("5 inputs required: InputData, Weights, N_Bins, Min_Value, Max_Value.");
   }

   mxClassID cid_inputdata = mxGetClassID(prhs[0]);
   int       dim_inputdata = mxGetNumberOfDimensions(prhs[0]); 

   mxClassID cid_weights   = mxGetClassID(prhs[1]);
   int       dim_weights   = mxGetNumberOfDimensions(prhs[1]); 

   int N;                  // Length of the input data array. 
   int N_classes;          // Number of classes: weights express the probability that each input value belongs to any of the N_classes. 
   int status = 1;         // Return 0 if everython ok, 1 if errors. 
   int N_bins;             
   double min_value;          
   double max_value;          

   /* The inputs must be noncomplex single or double floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Activity' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Cameras' must be noncomplex single or double.");

   if (dim_inputdata != 2 || mxGetDimensions(prhs[0])[1]!= 1)
       mexErrMsgTxt("'input_data' must be an array of size [Nx1]. ");

   if (dim_weights != 2)
       mexErrMsgTxt("'weights' must be an array of size [NxN_classes]. ");

   if (mxGetDimensions(prhs[0])[0] != mxGetDimensions(prhs[1])[0])
       mexErrMsgTxt("'input_data' and 'weights' must have the same number of elements: 'input_data':[Nx1]; 'weights':[NxN_classes].");
 
   /*  */
   N = mxGetDimensions(prhs[0])[0]; 
   N_classes = mxGetDimensions(prhs[1])[1]; 
   N_bins = (int) (mxGetScalar(prhs[2]));
   min_value = (double) (mxGetScalar(prhs[3])); 
   max_value = (double) (mxGetScalar(prhs[4])); 

   float *inputdata_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       inputdata_ptr = (float *) (mxGetData(prhs[0]));
   else
   {
       double *inputdata_ptr_d = (double *) (mxGetData(prhs[0]));
       inputdata_ptr = (float*) malloc(N*sizeof(float)); 
       for (int i=0; i<N;i++)
           inputdata_ptr[i] = inputdata_ptr_d[i];  
   }

   float *weights_ptr;
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       weights_ptr = (float *) (mxGetData(prhs[1]));
   else
   {
       double *weights_ptr_d = (double *) (mxGetData(prhs[1]));
       weights_ptr = (float*) malloc(N*N_classes*sizeof(float)); 
       for (int i=0; i<N*N_classes;i++)
           weights_ptr[i] = weights_ptr_d[i];  
   }

   /* Allocate array for weighted histogram */   
   mwSize mw_hist_size[2];
   mw_hist_size[0] = (mwSize)N_bins;
   mw_hist_size[1] = (mwSize)N_classes;

   plhs[0] =  mxCreateNumericArray(2, mw_hist_size, mxSINGLE_CLASS, mxREAL); 
   float *histogram_ptr = (float *)(mxGetData(plhs[0]));

   /* Calculate weighted histogram */
   status = et_array_histogram_weighted(inputdata_ptr, weights_ptr, histogram_ptr, N, N_classes, N_bins, min_value, max_value); 

   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(inputdata_ptr); 
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(weights_ptr); 
}


