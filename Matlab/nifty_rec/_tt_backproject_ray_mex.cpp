/*
 *  _tt_backproject_ray_mex.cpp
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

#define DEFAULT_USEGPU 0
#define DEFAULT_INTERP 0

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments. 
   if (!(nrhs==8 || nrhs==9 || nrhs==10) ){
      mexErrMsgTxt("8,9 or 10 arguments required: projections, detector_size, detector_translation, detector_rotation, source_position, volume_voxels, volume_size, t_step, [interpolation], [USE_GPU]");
   }

   int status = 1;
   float *projections;
   u_int_3 volume_voxels; 
   u_int   n_projections; 
   u_int_2 detector_pixels;
   float_2 *detector_size; 
   float_3 volume_size;
   float_3 *detector_translation;
   float_3 *detector_rotation;  
   float_3 *source_pos; 
   int enable_gpu=DEFAULT_USEGPU;
   int interpolation=DEFAULT_INTERP;

   //check consistency of arguments
   if ( mxGetClassID(prhs[0]) != mxSINGLE_CLASS ) mexErrMsgTxt("'projections' must be noncomplex single.");
   if ( mxGetClassID(prhs[1]) != mxDOUBLE_CLASS ) mexErrMsgTxt("'detector_size' must be noncomplex  double.");
   if ( mxGetClassID(prhs[2]) != mxDOUBLE_CLASS ) mexErrMsgTxt("'detector_translation' must be noncomplex  double.");
   if ( mxGetClassID(prhs[3]) != mxDOUBLE_CLASS ) mexErrMsgTxt("'detector_rotation' must be noncomplex  double.");
   if ( mxGetClassID(prhs[4]) != mxDOUBLE_CLASS ) mexErrMsgTxt("'source_position' must be noncomplex  double.");
   if ( mxGetClassID(prhs[5]) != mxINT32_CLASS )  mexErrMsgTxt("'volume_voxels' must be noncomplex  int32.");
   if ( mxGetClassID(prhs[6]) != mxDOUBLE_CLASS ) mexErrMsgTxt("'volume_size' must be noncomplex  double.");
   if ( mxGetClassID(prhs[7]) != mxDOUBLE_CLASS ) mexErrMsgTxt("'t_step' must be noncomplex  double.");

   if ( (mxGetNumberOfDimensions(prhs[0]) != 2) && (mxGetNumberOfDimensions(prhs[0]) != 3))   mexErrMsgTxt("'projections' must be a 2D or 3D matrix.");  
   if ( mxGetNumberOfDimensions(prhs[1]) != 2 )   mexErrMsgTxt("'detector_size' must be a 2D matrix.");  
   if ( mxGetNumberOfDimensions(prhs[2]) != 2 )   mexErrMsgTxt("'detector_translation' must be a 2D matrix."); 
   if ( mxGetNumberOfDimensions(prhs[3]) != 2 )   mexErrMsgTxt("'detector_rotation' must be a 2D matrix."); 
   if ( mxGetNumberOfDimensions(prhs[4]) != 2 )   mexErrMsgTxt("'source_position' must be a 2D matrix."); 
   if ( mxGetNumberOfDimensions(prhs[5]) != 2 )   mexErrMsgTxt("'volume_voxels' must be a 2D matrix."); 
   if ( mxGetNumberOfDimensions(prhs[6]) != 2 )   mexErrMsgTxt("'volume_size' must be a 2D matrix."); 
//   if ( mxGetNumberOfDimensions(prhs[7]) != 1 )   mexErrMsgTxt("'t_step' must be a scalar."); 

   if ( mxGetDimensions(prhs[2])[0] != mxGetDimensions(prhs[3])[0] || mxGetDimensions(prhs[3])[0] != mxGetDimensions(prhs[4])[0])  
       mexErrMsgTxt("size of parameters 3,4,5 should be the same: [N_projections,3]"); 
   if ( mxGetDimensions(prhs[2])[1] != mxGetDimensions(prhs[3])[1] || mxGetDimensions(prhs[3])[1] != mxGetDimensions(prhs[4])[1] ||  mxGetDimensions(prhs[4])[1] != 3)  
       mexErrMsgTxt("size of parameters 3,4,5 should be the same: [N_projections,3]"); 
   if ( mxGetDimensions(prhs[1])[0] != mxGetDimensions(prhs[2])[0] )
       mexErrMsgTxt("size of parameter 2 should be the same: [N_projections,2]");  
   if ( mxGetDimensions(prhs[1])[1] != 2 ) 
       mexErrMsgTxt("size of parameter 2 should be the same: [N_projections,2]");  
   if ( (mxGetDimensions(prhs[5])[1] != 3) | (mxGetDimensions(prhs[5])[0] != 1)) 
       mexErrMsgTxt("size of parameter 6 should be the same: [1,3]");  
   if ( (mxGetDimensions(prhs[6])[1] != 3) | (mxGetDimensions(prhs[6])[0] != 1)) 
       mexErrMsgTxt("size of parameter 7 should be the same: [1,3]");  

   n_projections = mxGetDimensions(prhs[2])[0];

   detector_pixels.w = mxGetDimensions(prhs[0])[0];
   detector_pixels.h = mxGetDimensions(prhs[0])[1];
   double* dsize_double = (double *) mxGetData(prhs[1]);
   double* dtranslation_double = (double *) mxGetData(prhs[2]);
   double* drotation_double = (double *) mxGetData(prhs[3]);
   double* sp_double   = (double *) mxGetData(prhs[4]);
   volume_voxels.x = ((int *) (mxGetData(prhs[5])))[0];
   volume_voxels.y = ((int *) (mxGetData(prhs[5])))[1];
   volume_voxels.z = ((int *) (mxGetData(prhs[5])))[2];
   volume_size.x = ((double *)mxGetData(prhs[6]))[0];
   volume_size.y = ((double *)mxGetData(prhs[6]))[1];
   volume_size.z = ((double *)mxGetData(prhs[6]))[2];
   float t_step = (mxGetScalar(prhs[7])); 
   if (nrhs>=9)
       interpolation = (int) (mxGetScalar(prhs[8]));
   if (nrhs>=10)
       enable_gpu = (int) (mxGetScalar(prhs[9]));


   //convert to CUDA compatible types
//   projections = (unsigned short*) malloc(detector_pixels.w*detector_pixels.h*n_projections*sizeof(float));
//   for(int i=0; i<detector_pixels.w*detector_pixels.h*n_projections; i++)
//       projections[i] = proj_double[i];
   projections = (float *) mxGetData(prhs[0]);

   detector_size        = (float_2 *) malloc(n_projections*sizeof(float_2));
   detector_translation = (float_3 *) malloc(n_projections*sizeof(float_3));
   detector_rotation    = (float_3 *) malloc(n_projections*sizeof(float_3));
   source_pos           = (float_3 *) malloc(n_projections*sizeof(float_3));
   for(int i=0;i<n_projections;i++)
   {
       source_pos[i].x     = sp_double[i];   source_pos[i].y     = sp_double[i+n_projections];   source_pos[i].z     = sp_double[i+2*n_projections];
       detector_size[i].x = dsize_double[i]; detector_size[i].y = dsize_double[i+n_projections]; 
       detector_translation[i].x = dtranslation_double[i]; detector_translation[i].y = dtranslation_double[i+n_projections]; detector_translation[i].z = dtranslation_double[i+2*n_projections];
       detector_rotation[i].x = drotation_double[i]; detector_rotation[i].y = drotation_double[i+n_projections]; detector_rotation[i].z = drotation_double[i+2*n_projections];
   }

   //allocate memory for backprojection
//   float *out_backprojection_float = (float*) malloc(volume_voxels.x*volume_voxels.y*volume_voxels.z*sizeof(float));
   mwSize mw_out_size[3];
   mw_out_size[0] = (mwSize)volume_voxels.x;
   mw_out_size[1] = (mwSize)volume_voxels.y;
   mw_out_size[2] = (mwSize)volume_voxels.z;
//   plhs[0] =  mxCreateNumericArray(3, mw_out_size, mxDOUBLE_CLASS, mxREAL);
//   double *out_double = (double *)(mxGetData(plhs[0]));   
   plhs[0] =  mxCreateNumericArray(3, mw_out_size, mxSINGLE_CLASS, mxREAL);
   float *out_single = (float *)(mxGetData(plhs[0]));   

fprintf(stderr,"---MEX> %f %f %f - %d %d %d - %f %f \n",volume_size.x,volume_size.y,volume_size.z,volume_voxels.x,volume_voxels.y,volume_voxels.z,detector_size[0].x,detector_size[0].y);

   // Backproject   
   status = tt_array_backproject_ray(projections, detector_pixels, n_projections, out_single, detector_size, detector_translation, detector_rotation, source_pos, volume_voxels, volume_size, t_step, interpolation, enable_gpu);


   // Convert to double 
//   for (int i=0; i<n_projections*detector_pixels.w*detector_pixels.h; i++)
//       out_double[i] = out_projections_float[i];

   // Shutdown 
//   free(attenuation);
   free(source_pos);
   free(detector_size);
   free(detector_translation);
   free(detector_rotation);
//   free(out_projections_float);

   // Return 
   if (status != 0)
   	mexErrMsgTxt("Error while performing Affine Transform.");
   return;
}


