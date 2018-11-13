/*
 *  _et_reset_gpu_mex.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include <mex.h>
#include "_et_array_interface.h"

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   
   /* Reset GPU */
   int status = et_array_reset_gpu();

   mwSize mw_status_size[1];
   mw_status_size[0] = (mwSize)1;
   plhs[0] =  mxCreateNumericArray(1, mw_status_size, mxINT32_CLASS, mxREAL);
   int *status_ptr = (int *)(mxGetData(plhs[0]));
   status_ptr[0]=status;

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("NiftyRec: Error while resetting the GPU.");
   return;
}


