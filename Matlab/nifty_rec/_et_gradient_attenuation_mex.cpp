/*
 *  _et_gradient_attenuation_mex.cpp
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
#include <math.h>
#include <cmath>


/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (!(nrhs==5 || nrhs==6 || nrhs==7 || nrhs==8)){
      mexErrMsgTxt("5,6,7 or 8 inputs required: Activity, Sinogram, Cameras, Attenuation, [PointSpreadFunction], [EnableGpu], [Background], [BackgroundAttenuation], [TruncateNegativeValues]");
   }

   mxClassID cid_activity = mxGetClassID(prhs[0]);
   int       dim_activity = mxGetNumberOfDimensions(prhs[0]); 

   mxClassID cid_sino = mxGetClassID(prhs[1]);
   int       dim_sino = mxGetNumberOfDimensions(prhs[1]); 
   
   mxClassID cid_cameras = mxGetClassID(prhs[2]);
   int       dim_cameras = mxGetNumberOfDimensions(prhs[2]);

   mxClassID cid_attenuation = mxGetClassID(prhs[3]);
   int       dim_attenuation = mxGetNumberOfDimensions(prhs[3]);

   mxClassID cid_psf = mxGetClassID(prhs[4]);
   int       dim_psf = mxGetNumberOfDimensions(prhs[4]);   

   int activity_size[3]; // Size of input activity matrix
   int sino_size[3];     // Size of input sinogram matrix
   int psf_size[3];      // Size of input psf matrix
   int gradient_size[3];     // Size of output backprojection matrix
   int cameras_size[2];  // Size of cameras matrix (can be (nx3) or (1xn))
   int enable_gpu = 0;      // Flag that enables(1)/disables(0) GPU acceleration
   float background = 0; // Background (for rotation in the backprojection)
   float background_attenuation=0; // Attenuation background (when rotating attenuation)
   int no_psf = 0;       // Flag for psf: if 1 it means that no PSF was specified
   int truncate_negative_values = 1; // Set this flag to 1 in order to truncate the results to 0 if negative. Set it to 0 to disable truncation. 

   int status = 1;       // Return status: 0 if succesful

   /* The inputs must be noncomplex single floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Activity' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Sinogram' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Cameras' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Attenuation' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[4]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Psf' must be noncomplex single or double.");

   /* Check if size of cameras matrix is correct */     
   cameras_size[0] = mxGetDimensions(prhs[2])[0];
   cameras_size[1] = mxGetDimensions(prhs[2])[1];
   if (!(cameras_size[1] == 3 || cameras_size[1] == 1))
      mexErrMsgTxt("Cameras must be of size [n_cameras x 1] or [n_cameras x 3]");     
   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: Backprojection");
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

   /* Check consistency of input (and create size of sinogram for return */
   switch(dim_sino)
      {
      /* Check consistency of input if 2D (and create size of sinogram for return) */
      case 2:
           if (!no_psf)
               if (dim_psf != 2)
                   mexErrMsgTxt("Dimension of Point Spread Function matrix must match the simension of Sinogram (ddpsf).");        
           sino_size[0] = mxGetDimensions(prhs[1])[0];
           sino_size[1] = mxGetDimensions(prhs[1])[1];
           sino_size[2] = 1;
           if (sino_size[0]<2)
               mexErrMsgTxt("Size of Activity matrix must be greater then 2");
           if (sino_size[1] != cameras_size[0])
               mexErrMsgTxt("Number of cameras must match Sinogram size");    
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[4])[0];
               psf_size[1] = mxGetDimensions(prhs[4])[1];
               psf_size[2] = 1;
               if (psf_size[0]%2!=1 || psf_size[1]!=sino_size[0])
                   mexErrMsgTxt("Point Spread Function must be of size hxN for Backprojection of size Nxn; h odd.");
               }
           gradient_size[0] = sino_size[0];
           gradient_size[1] = sino_size[0];
           gradient_size[2] = 1;
           if(mxGetDimensions(prhs[1])[0] != mxGetDimensions(prhs[3])[0] || mxGetDimensions(prhs[1])[1] != mxGetDimensions(prhs[3])[1])
               mexErrMsgTxt("Attenuation must be of the same size of Activity");
           break;
      /* Check consistency of input if 3D (and create size of sinogram for return */
      case 3:
           if (!no_psf)
               if (dim_psf != 3)
                   mexErrMsgTxt("Dimension of Point Spread Function matrix must match the dimension of Sinogram (ddpsf).");           
           sino_size[0] = mxGetDimensions(prhs[1])[0];
           sino_size[1] = mxGetDimensions(prhs[1])[1];
           sino_size[2] = mxGetDimensions(prhs[1])[2];
//           if (sino_size[0]<2 || sino_size[1]<2)
//               mexErrMsgTxt("Size of Activity matrix must be greater then 2");           
           if (sino_size[2] != cameras_size[0])
                mexErrMsgTxt("Number of cameras must match Sinogram size");           
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[4])[0];
               psf_size[1] = mxGetDimensions(prhs[4])[1];
               psf_size[2] = mxGetDimensions(prhs[4])[2];
               if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=sino_size[0])
                    mexErrMsgTxt("Point Spread Function must be of size hxkxN for Activity of size NxmxN; h,k odd.");
               }
           gradient_size[0] = sino_size[0];
           gradient_size[1] = sino_size[1];
           gradient_size[2] = sino_size[0];
           activity_size[0] = mxGetDimensions(prhs[0])[0];
           activity_size[1] = mxGetDimensions(prhs[0])[1];
           activity_size[2] = mxGetDimensions(prhs[0])[2]; 
           if ( activity_size[0]!=gradient_size[0] || activity_size[1]!=gradient_size[1] || activity_size[2]!=gradient_size[2])
               mexErrMsgTxt("Activity must be of size [NxMxN] for sinogram of size [NxMxNcameras]");
           if ((activity_size[0]!=activity_size[2]))
               mexErrMsgTxt("3D Activity must be of size [NxMxN]");
           if(gradient_size[0] != mxGetDimensions(prhs[3])[0] || gradient_size[1] != mxGetDimensions(prhs[3])[1] || gradient_size[2] != mxGetDimensions(prhs[3])[2])
               mexErrMsgTxt("Attenuation must be of the same size of Activity");
           break;        
      default:
           mexErrMsgTxt("Activity must be either 2D or 3D.");
           break;
      }

   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>=6)
       enable_gpu = (int) (mxGetScalar(prhs[5]));

   /* Check if activity is multiple of ET_BLOCK_SIZE */
   if (enable_gpu) {
       if (!et_is_block_multiple(gradient_size[0]) || !et_is_block_multiple(gradient_size[1])) {
           char msg[100];
           sprintf(msg,"With GPU enabled, size of activity must be a multiple of %d",et_get_block_size());
           mexErrMsgTxt(msg);
           }
       }

   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=7)
       background = (mxGetScalar(prhs[6]));

   /* Check if BackgroundAttenuation is specified */
   background_attenuation = 0.0f;
   if (nrhs>=8)
       background_attenuation = (mxGetScalar(prhs[7]));

   /* Check if 'truncate negative values' flag has been specified */ 
   if (nrhs>=9)
       truncate_negative_values = (mxGetScalar(prhs[8])); 

   /* Extract pointers to input matrices */
   float *sinogram_ptr;
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       sinogram_ptr = (float *) (mxGetData(prhs[1]));
   else
   {
       double *sinogram_ptr_d = (double *) (mxGetData(prhs[1]));
       sinogram_ptr = (float*) malloc(sino_size[0]*sino_size[1]*sino_size[2]*sizeof(float));
       for (int i=0; i<sino_size[0]*sino_size[1]*sino_size[2];i++)
           sinogram_ptr[i] = sinogram_ptr_d[i];
   }

   float *cameras_ptr;
   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       cameras_ptr = (float *) (mxGetData(prhs[2]));
   else
   {
       double *cameras_ptr_d = (double *) (mxGetData(prhs[2]));
       cameras_ptr = (float*) malloc(cameras_size[0]*cameras_size[1]*sizeof(float));
       for (int i=0; i<cameras_size[0]*cameras_size[1];i++)
           cameras_ptr[i] = cameras_ptr_d[i];  
   }

   int attenuation_size[3];
   float *attenuation_ptr=NULL;
   attenuation_size[0] = gradient_size[0];
   attenuation_size[1] = gradient_size[1];
   attenuation_size[2] = gradient_size[2];
   if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
       attenuation_ptr = (float *) (mxGetData(prhs[3]));
   else
       {
       double *attenuation_ptr_d = (double *) (mxGetData(prhs[3]));
       attenuation_ptr = (float *) malloc(attenuation_size[0]*attenuation_size[1]*attenuation_size[2]*sizeof(float));
       for (int i=0; i<attenuation_size[0]*attenuation_size[1]*attenuation_size[2];i++)
           attenuation_ptr[i] = attenuation_ptr_d[i];  
       }

   float *activity_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       activity_ptr = (float *) (mxGetData(prhs[0]));
   else
   {
       double *activity_ptr_d = (double *) (mxGetData(prhs[0]));
       activity_ptr = (float*) malloc(activity_size[0]*activity_size[1]*activity_size[2]*sizeof(float));
       for (int i=0; i<activity_size[0]*activity_size[1]*activity_size[2];i++)
           activity_ptr[i] = activity_ptr_d[i];
   }

   float *psf_ptr=NULL;
   if (!no_psf)
       {
       if (mxGetClassID(prhs[4]) == mxSINGLE_CLASS)
           psf_ptr = (float *) (mxGetData(prhs[4]));
       else
           {
           double *psf_ptr_d = (double *) (mxGetData(prhs[4]));
           psf_ptr = (float *) malloc(psf_size[0]*psf_size[1]*psf_size[2]*sizeof(float));
           for (int i=0; i<psf_size[0]*psf_size[1]*psf_size[2];i++)
               psf_ptr[i] = psf_ptr_d[i];  
           }
       }
              
   /* Allocate gradient matrix */   
   int dim_gradient;   
   mwSize mw_gradient_size[3];
   
   dim_gradient = dim_sino;
   mw_gradient_size[0] = (mwSize)gradient_size[0];
   mw_gradient_size[1] = (mwSize)gradient_size[1]; 
   mw_gradient_size[2] = (mwSize)gradient_size[2]; 

   plhs[0] =  mxCreateNumericArray(dim_gradient, mw_gradient_size, mxSINGLE_CLASS, mxREAL);
   float *gradient_ptr = (float *)(mxGetData(plhs[0]));

   /* Calculate gradient */
   status = et_array_gradient_attenuation(sinogram_ptr, sino_size, activity_ptr, activity_size, gradient_ptr, gradient_size, cameras_ptr, cameras_size, psf_ptr, psf_size, attenuation_ptr, attenuation_size, background, background_attenuation, enable_gpu, truncate_negative_values);

   /* Shutdown */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(activity_ptr); 
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(sinogram_ptr);
   if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) free(cameras_ptr);
   if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) free(attenuation_ptr); 
   if ((!no_psf) && (mxGetClassID(prhs[4]) != mxSINGLE_CLASS)) free(psf_ptr);
   
   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while calculating gradient.");
   return;
}


