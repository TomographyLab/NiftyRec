
/*
 *  _et_list_gpus_mex.cpp
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


/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs>0){
      mexErrMsgTxt("No arguments expected");
   }

   int status = 1;

   status = et_array_isinstalled();
   if (status != 0)
        {
   	mexErrMsgTxt("NiftyRec is not working, please reinstall.");
        }

   /* Found some GPU devices: return matrix with GPUs info */   
   mwSize mw_size[1];
   mw_size[0] = (mwSize)1;
   plhs[0] =  mxCreateNumericArray(1, mw_size, mxINT32_CLASS, mxREAL);
   int *return_ptr = (int *)(mxGetData(plhs[0]));
   if (status==0) 
       return_ptr[0] = 1; 
   else
       return_ptr[0] = 0;
   return;
}


