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
   int gpu_count;
   int *gpus_info_array = (int*) malloc(MAX_DEVICES*SIZE_OF_INFO*sizeof(int));

   status = et_array_list_gpus(&gpu_count, gpus_info_array);
   if (status != 0)
        {
        free(gpus_info_array);
   	mexErrMsgTxt("Error while querying GPU info.");
        }

   if (gpu_count == 0)
       {
       /* No GPUs: return 0 */  
       fprintf(stderr,"\net_list_gpus: No GPU support");
       mwSize mw_info_size[1];
       mw_info_size[0] = (mwSize)1;
       plhs[0] =  mxCreateNumericArray(1, mw_info_size, mxINT32_CLASS, mxREAL);
       int *info_ptr = (int *)(mxGetData(plhs[0]));
       info_ptr[0]=0;
       }
   else
       {
       /* Found some GPU devices: return matrix with GPUs info */   
       mwSize mw_info_size[2];
       mw_info_size[0] = (mwSize)gpu_count;
       mw_info_size[1] = (mwSize)SIZE_OF_INFO; 
       plhs[0] =  mxCreateNumericArray(2, mw_info_size, mxINT32_CLASS, mxREAL);
       int *info_ptr = (int *)(mxGetData(plhs[0]));

       /* Copy to Matlab matrix (and reorder for Matlab) */
       for (int gpu=0; gpu<gpu_count; gpu++)
           {
           for (int field=0; field<SIZE_OF_INFO; field++)
               info_ptr[gpu_count*field+gpu] = gpus_info_array[SIZE_OF_INFO*gpu+field];
           }
       }
   free(gpus_info_array);
   return;
}


