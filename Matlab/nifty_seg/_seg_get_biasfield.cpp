/*
 *  _seg_get_biasfield.cpp
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

   /* Instanciate return matrix for biasfield data */
   int size_x=0;
   int size_y=0;
   int size_z=0;

   status = seg_array_get_biasfield_size(&size_x, &size_y, &size_z);
   fprintf(stderr,"Bias Field size: %d %d %d \n",size_x, size_y, size_z);

   mwSize mw_out_size[3]; //only works with 3D images, FIXME
   mw_out_size[0] = (mwSize)size_x;
   mw_out_size[1] = (mwSize)size_y;
   mw_out_size[2] = (mwSize)size_z;

   plhs[0] =  mxCreateNumericArray(3, mw_out_size, mxSINGLE_CLASS, mxREAL); //not templated, only works with float FIXME
   PrecisionTYPE* out_ptr = (PrecisionTYPE *)(mxGetData(plhs[0]));

   status = seg_array_get_biasfield(out_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while retrieving bias field corrected image.");
   return;
}


