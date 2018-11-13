/*
 *  _seg_get_mean_variance.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include <mex.h>
#include "_seg_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int status;

   /* Check for proper number of arguments. */
   if (!(nrhs==0)){
      mexErrMsgTxt("No input parameters expected");
   }

   /* Instanciate return matrix for segmentation data */
   int size_x=0;
   int size_y=0;
   int size_z=0;
   int n_classes=0; 
   int n_images=0;

   status = seg_array_get_segmentation_size(&size_x, &size_y, &size_z, &n_classes);
   if (status != 0)
   {
   	mexErrMsgTxt("Error while retrieving priors.");
       return;
   }

   status = seg_array_get_image_size(&size_x, &size_y, &size_z, &n_images);
   if (status != 0)
   {
   	mexErrMsgTxt("Error while retrieving priors.");
       return;
   }

   mwSize mw_mean_size[2];
   mwSize mw_var_size[3];

   if (n_images == 1)
       {
       mw_mean_size[0] = (mwSize)1;
       mw_mean_size[1] = (mwSize)n_classes;
       plhs[0] =  mxCreateNumericArray(2, mw_mean_size, mxDOUBLE_CLASS, mxREAL); 
       mw_var_size[0] = (mwSize)1;
       mw_var_size[1] = (mwSize)n_classes;
       plhs[1] =  mxCreateNumericArray(2, mw_var_size, mxDOUBLE_CLASS, mxREAL); 
       }
   else
       {
       mw_mean_size[0] = (mwSize)n_images;
       mw_mean_size[1] = (mwSize)n_classes;
       plhs[0] =  mxCreateNumericArray(2, mw_mean_size, mxDOUBLE_CLASS, mxREAL); 
       mw_var_size[0] = (mwSize)n_images;
       mw_var_size[1] = (mwSize)n_images;
       mw_var_size[2] = (mwSize)n_classes;
       plhs[1] =  mxCreateNumericArray(3, mw_var_size, mxDOUBLE_CLASS, mxREAL); 
       }

   double* mean_ptr = (double *)(mxGetData(plhs[0]));
   double* var_ptr = (double *)(mxGetData(plhs[1]));

   status = seg_array_get_mean_variance(mean_ptr, var_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while retrieving the segmentation.");
   return;
}


