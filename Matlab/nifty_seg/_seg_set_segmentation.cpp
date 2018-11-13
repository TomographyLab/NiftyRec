/*
 *  _seg_set_segmentation.cpp
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
   if (!(nrhs==1))
      mexErrMsgTxt("1 input required: probabilistic segmentation");

   /* Check for correct size and type */
   const int dim_in = mxGetNumberOfDimensions(prhs[0]); 
   if (dim_in!=4)
      mexErrMsgTxt("probabilistic segmentation must be of size (size_x, size_y, size_z, n_classes)");

   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
      mexErrMsgTxt("probabilistic segmentation must be of type SINGLE. Use  single()");

   int size_x=0, size_x_in=0;
   int size_y=0, size_y_in=0;
   int size_z=0, size_z_in=0;
   int n_classes=0, n_classes_in=0;

   status = seg_array_get_segmentation_size(&size_x, &size_y, &size_z, &n_classes);
   fprintf(stderr,"Segmentation size: %d %d %d %d \n",size_x, size_y, size_z, n_classes);

   size_x_in = mxGetDimensions(prhs[0])[0];
   size_y_in = mxGetDimensions(prhs[0])[1];
   size_z_in = mxGetDimensions(prhs[0])[2];
   n_classes_in = mxGetDimensions(prhs[0])[3];

   if (size_x!=size_x_in || size_y!=size_y_in || size_z!=size_z_in || n_classes!=n_classes_in)
       mexErrMsgTxt("size of probabilistic segmentation must match the size of the input data and the number of classes");

   /* Call NiftySeg function */
   float *data_ptr = (float*) mxGetData(prhs[0]); //not templated, FIXME
   status = seg_array_set_segmentation(data_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("NiftySeg: error while setting the probabilistic segmentation.");
   return;
}


