/*
 *  _tt_project_perspective_mex.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <mex.h>

#include "_tt_array_interface.h"

#include <limits>
#include <string.h>
#include <math.h>
#include <cmath>

#define eps 0.00001f

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (!(nrhs==5 || nrhs==5 || nrhs==6 || nrhs==7)){
      mexErrMsgTxt("4, 5, 6 or 7 inputs required: Attenuation, ImageOrigin, DetectorOrigin, DetectorShape, PointSpreadFunction, EnableGPU, Background");
   }

   mxClassID cid_psf;
   int       dim_psf          = 0;

   int attenuation_size[3];  //size of input attenuation.
   int psf_size[3];          //size of input psf.
   int projection_size[3];   //size of output projection.
   int detector_origin_size[2];  
   int image_origin_size[2]; 
   int detector_shape_size[2]; 
   float background;         // Background value, defaults to 0. (It is used to fill empty areas when images are rotated)
   int enable_gpu;           // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int no_psf = 0;           // This flag goes high if psf input parameter from the Matlab function is a scalar -> no psf

   int N;                    // Size of attenuation: [N,N,m].
   int m;                    // Size of attenuation: [N,N,m].
   int status = 1;           // Return 0 if everything ok, 1 if errors.

   /* The inputs must be noncomplex single floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Attenuation' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'ImageOrigin' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'DetectorOrigin' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'DetectorShape' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'PointSpreadFunction' must be noncomplex single or double.");

   /* Check if size of detectors matrix is correct */
   image_origin_size[0] = mxGetDimensions(prhs[1])[0];
   image_origin_size[1] = mxGetDimensions(prhs[1])[1];
   detector_origin_size[0] = mxGetDimensions(prhs[2])[0];
   detector_origin_size[1] = mxGetDimensions(prhs[2])[1];
   detector_shape_size[0] = mxGetDimensions(prhs[3])[0];
   detector_shape_size[1] = mxGetDimensions(prhs[3])[1];
   if (!(image_origin_size[1] == 3))
      mexErrMsgTxt("ImageOrigin must be of size [n_projections x 3]"); 
   if (!(detector_origin_size[1] == 3))
      mexErrMsgTxt("DetectorOrigin must be of size [n_projections x 3]"); 
   if (!(detector_shape_size[1] == 2))
      mexErrMsgTxt("DetectorShape must be of size [n_projections x 2]");
   if (! (image_origin_size[0]==detector_origin_size[0] && detector_origin_size[0]==detector_shape_size[0]) ) 
      mexErrMsgTxt("DetectorOrigin, ImageOrigin and DetectorShape must be of size: n_projections");
   int n_projections = image_origin_size[0];

   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: Projection");
   }     
   /* Check is psf is a scalar, in that case do not apply any psf */
   if (nrhs<5)  //no psf parameter
       no_psf = 1;
   else
       {
       cid_psf = mxGetClassID(prhs[4]);
       dim_psf = mxGetNumberOfDimensions(prhs[4]);  
       int all_one=1;
       for (int i=0; i<dim_psf; i++)
           if (mxGetDimensions(prhs[4])[i]!=1)
              all_one = 0;
       if (all_one)
           no_psf = 1;  //psf parameter is a scalar
       }
   if (no_psf == 1)
       {
       psf_size[0] = 0;
       psf_size[1] = 0;
       psf_size[2] = 0;
       }
   /* Check consistency of input (and create size of attenuation for return */
   if (!no_psf)
          if (dim_psf != 3 && !no_psf)
                  mexErrMsgTxt("Size of Point Spread Function matrix must match the size of Attenuation (ddpsf).");           
   attenuation_size[0] = mxGetDimensions(prhs[0])[0];
   attenuation_size[1] = mxGetDimensions(prhs[0])[1];
   attenuation_size[2] = mxGetDimensions(prhs[0])[2];
   if (!no_psf)
         {
         psf_size[0] = mxGetDimensions(prhs[4])[0];
         psf_size[1] = mxGetDimensions(prhs[4])[1];
         psf_size[2] = mxGetDimensions(prhs[4])[2];
         if ((psf_size[0]%2!=1) || ((psf_size[1]%2!=1)) || psf_size[2]!=attenuation_size[2])
             mexErrMsgTxt("Point Spread Function must be of size hxkxm for attenuation of size NxNxm; h,k odd");
         }
   N = attenuation_size[0];
   m = attenuation_size[2];
   if (attenuation_size[0]<2 || attenuation_size[1]<2 || attenuation_size[2]<2)
         mexErrMsgTxt("Size of attenuation matrix must be greater then 2");           
   if (attenuation_size[1]!=attenuation_size[0])
         mexErrMsgTxt("2D attenuation must be square (3D attenuation can be of size NxNxm)");
   projection_size[0] = attenuation_size[0];
   projection_size[1] = attenuation_size[2];
   projection_size[2] = detector_origin_size[0];      


   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>=4)
       enable_gpu = (int) (mxGetScalar(prhs[5]));
  
   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=5)
       background = (mxGetScalar(prhs[6])); 

   /* Extract pointers to input matrices (and eventually convert double to float)*/
   float *temp_vector = (float*) malloc(n_projections*3*sizeof(float));

   float *attenuation_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       attenuation_ptr = (float *) (mxGetData(prhs[0]));
   else
   {
       double *attenuation_ptr_d = (double *) (mxGetData(prhs[0]));
       attenuation_ptr = (float*) malloc(attenuation_size[0]*attenuation_size[1]*attenuation_size[2]*sizeof(float));
       for (int i=0; i<attenuation_size[0]*attenuation_size[1]*attenuation_size[2];i++)
           attenuation_ptr[i] = attenuation_ptr_d[i];
   }

   float *image_origin_ptr;
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       image_origin_ptr = (float *) (mxGetData(prhs[1]));
   else
   {
       double *image_origin_ptr_d = (double *) (mxGetData(prhs[1]));
       image_origin_ptr = (float*) malloc(image_origin_size[0]*image_origin_size[1]*sizeof(float));
       for (int i=0; i<image_origin_size[0]*image_origin_size[1];i++)
           image_origin_ptr[i] = image_origin_ptr_d[i];  
   }
   //reorder for compatibility with c-array interface
   for(int i=0; i<n_projections; i++)
       {
       temp_vector[3*i] = image_origin_ptr[i];
       temp_vector[3*i+1] = image_origin_ptr[i+1*n_projections];
       temp_vector[3*i+2] = image_origin_ptr[i+2*n_projections];
       }
   for(int i=0; i<3*n_projections; i++)
       image_origin_ptr[i]=temp_vector[i];


   float *detector_origin_ptr;
   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       detector_origin_ptr = (float *) (mxGetData(prhs[2]));
   else
   {
       double *detector_origin_ptr_d = (double *) (mxGetData(prhs[2]));
       detector_origin_ptr = (float*) malloc(detector_origin_size[0]*detector_origin_size[1]*sizeof(float));
       for (int i=0; i<detector_origin_size[0]*detector_origin_size[1];i++)
           detector_origin_ptr[i] = detector_origin_ptr_d[i];  
   }

   //reorder for compatibility with c-array interface
   for(int i=0; i<n_projections; i++)
       {
       temp_vector[3*i] = detector_origin_ptr[i];
       temp_vector[3*i+1] = detector_origin_ptr[i+1*n_projections];
       temp_vector[3*i+2] = detector_origin_ptr[i+2*n_projections];
       }
   for(int i=0; i<3*n_projections; i++)
       detector_origin_ptr[i]=temp_vector[i];


   float *detector_shape_ptr;
   if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
       detector_shape_ptr = (float *) (mxGetData(prhs[3]));
   else
   {
       double *detector_shape_ptr_d = (double *) (mxGetData(prhs[3]));
       detector_shape_ptr = (float*) malloc(detector_shape_size[0]*detector_shape_size[1]*sizeof(float));
       for (int i=0; i<detector_shape_size[0]*detector_shape_size[1];i++)
           detector_shape_ptr[i] = detector_shape_ptr_d[i];  
   }
   //reorder for compatibility with c-array interface
   for(int i=0; i<n_projections; i++)
       {
       temp_vector[2*i] = detector_shape_ptr[i];
       temp_vector[2*i+1] = detector_shape_ptr[i+1*n_projections];
       }
   for(int i=0; i<2*n_projections; i++)
       detector_shape_ptr[i]=temp_vector[i];


   float *psf_ptr;
   if (mxGetClassID(prhs[4]) == mxSINGLE_CLASS)
       psf_ptr = (float *) (mxGetData(prhs[4]));
   else
   {
       double *psf_ptr_d = (double *) (mxGetData(prhs[4]));
       psf_ptr = (float *) malloc(psf_size[0]*psf_size[1]*psf_size[2]*sizeof(float));
       for (int i=0; i<psf_size[0]*psf_size[1]*psf_size[2];i++)
           psf_ptr[i] = psf_ptr_d[i];  
   }

   /* Allocate projection matrix */     
   mwSize mw_projection_size[3];
   
   mw_projection_size[0] = (mwSize)projection_size[0];
   mw_projection_size[1] = (mwSize)projection_size[1]; 
   mw_projection_size[2] = (mwSize)projection_size[2]; 

   plhs[0] =  mxCreateNumericArray(3, mw_projection_size, mxSINGLE_CLASS, mxREAL);
   float *projection_ptr = (float *)(mxGetData(plhs[0]));

   /* Perform projection */
   status = tt_array_project_perspective(attenuation_ptr, attenuation_size, projection_ptr, projection_size, image_origin_ptr, detector_origin_ptr, detector_shape_ptr, psf_ptr, psf_size, background, enable_gpu);

   /* Shutdown */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(attenuation_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(image_origin_ptr);
   if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) free(detector_origin_ptr);
   if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) free(detector_shape_ptr);
   if (mxGetClassID(prhs[4]) != mxSINGLE_CLASS) free(psf_ptr);
   free(temp_vector);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("tt_project_mex: Error while performing projection.");
   return;
}


