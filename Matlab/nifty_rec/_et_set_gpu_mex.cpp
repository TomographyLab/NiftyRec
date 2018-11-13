/*
 *  _et_set_gpu_mex.cpp
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
   if (nrhs!=1){
      mexErrMsgTxt("1 scalar input is expected: DeviceID");
   }
   int id = (int) (mxGetScalar(prhs[0]));

   int status = 1;
   int gpu_count;
   int *gpus_info_array = (int*) malloc(MAX_DEVICES*SIZE_OF_INFO*sizeof(int));

   int matching_gpu = -1;

   status = et_array_list_gpus(&gpu_count, gpus_info_array);
   if (status != 0)
        {
        free(gpus_info_array);
   	mexErrMsgTxt("Error while querying GPU info.");
        }

   if (gpu_count == 0)
       {
       /* No GPUs: return 0 */  
       status = 1;
       }
   else
       {
       /* Found some GPU devices: check if the requested device ID is there */   
       for (int gpu=0; gpu<gpu_count; gpu++)
           {
           for (int field=0; field<SIZE_OF_INFO; field++)
               if (gpus_info_array[SIZE_OF_INFO*gpu+0] == id)
                   matching_gpu = gpu;
           }
       }
       if (matching_gpu>-1)
           fprintf(stderr,"\nHere");
           status = et_array_set_gpu(id);

   /* Could not set gpu device: return 0 */ 
   if (status==1)
       {
       fprintf(stderr,"\net_set_gpu: No GPU support");
       mwSize mw_info_size[1];
       mw_info_size[0] = (mwSize)1;
       plhs[0] =  mxCreateNumericArray(1, mw_info_size, mxINT32_CLASS, mxREAL);
       int *info_ptr = (int *)(mxGetData(plhs[0]));
       info_ptr[0]=0;
       }
   else
       {
       mwSize mw_info_size[2];
       mw_info_size[0] = 1;
       mw_info_size[1] = (mwSize)SIZE_OF_INFO; 
       plhs[0] = mxCreateNumericArray(2, mw_info_size, mxINT32_CLASS, mxREAL);
       int *info_ptr = (int *)(mxGetData(plhs[0]));
       /* Copy to Matlab matrix (and reorder for Matlab) */
       for (int field=0; field<SIZE_OF_INFO; field++)
           info_ptr[field] = gpus_info_array[SIZE_OF_INFO*matching_gpu+field];
       }
   free(gpus_info_array);
   return;
}


