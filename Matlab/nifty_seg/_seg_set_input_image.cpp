/*
 *  _seg_set_input_image.cpp
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
      mexErrMsgTxt("1 input required: InputImage");

   /* Check for correct size and type */
   const int dim_in = mxGetNumberOfDimensions(prhs[0]); 
   if (!(dim_in==3 || dim_in==4))
      mexErrMsgTxt("Input image must be of size (Nx,Ny,Nz) or (Nx,Ny,Nz,Nimages)");

   int size_x=1, size_x_in=1;
   int size_y=1, size_y_in=1;
   int size_z=1, size_z_in=1;
   int n_images=1, n_images_in=1;

   status = seg_array_get_image_size(&size_x, &size_y, &size_z, &n_images);

   size_x_in = mxGetDimensions(prhs[0])[0];
   size_y_in = mxGetDimensions(prhs[0])[1];
   size_z_in = mxGetDimensions(prhs[0])[2];
   if (dim_in>=4)
       n_images_in = mxGetDimensions(prhs[0])[3];
   else
       n_images_in = 1;

   if (size_x!=size_x_in || size_y!=size_y_in || size_z!=size_z_in || n_images!=n_images_in)
       mexErrMsgTxt("size of InputImage must match the size of the initialisation.");

   if (!(mxGetClassID(prhs[0]) == mxSINGLE_CLASS || mxGetClassID(prhs[0]) == mxDOUBLE_CLASS))
       mexErrMsgTxt("Input type must be single or double.");

   /* Extract input iamge data */
   float *data_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
      {
      data_ptr = (float*) mxGetData(prhs[0]);
      }
   else
      {
      double *data_ptr_d = (double*) mxGetData(prhs[0]); 
      data_ptr = (float*) malloc(size_x*size_y*size_z*n_images*sizeof(float));
      for (int i = 0 ; i<size_x*size_y*size_z*n_images; i++)
          data_ptr[i] = data_ptr_d[i];
      }

   /* Call NiftySeg function */
   status = seg_array_set_input_image(data_ptr);

   /* Return */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
       free(data_ptr);

   if (status != 0)
   	mexErrMsgTxt("NiftySeg: error while setting the probabilistic segmentation.");
   return;
}


