/*
 *  _seg_get_priors.cpp
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
//   if (!(nrhs==1)){
//      mexErrMsgTxt("1 input parameter: class_number");
//   }
   if (!(nrhs==0)){
      mexErrMsgTxt("No input parameters.");
   }

   /* Instanciate return matrix for segmentation data */
   int size_x=0;
   int size_y=0;
   int size_z=0;
   int n_images;
   int n_classes=0; 
//   int class_number=0;

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

//   class_number = (int) mxGetScalar(prhs[0]);

   mwSize mw_out_size[4];   // only 3D, FIXME
   mw_out_size[0] = (mwSize)size_x;
   mw_out_size[1] = (mwSize)size_y;
   mw_out_size[2] = (mwSize)size_z;
   mw_out_size[3] = (mwSize)n_classes;

   plhs[0] =  mxCreateNumericArray(4, mw_out_size, mxSINGLE_CLASS, mxREAL); //not templated, FIXME
   PrecisionTYPE* out_ptr = (PrecisionTYPE *)(mxGetData(plhs[0]));

   status = seg_array_get_priors(out_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while retrieving the priors, perhaps they have not been loaded.");
   return;
}


