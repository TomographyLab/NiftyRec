/*
 *  _et_backprojection_mex.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
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
   if (!(nrhs==4 || nrhs==5 || nrhs==6 || nrhs==7 || nrhs==8)){
      mexErrMsgTxt("2,3,4,5,6,7 or 8 inputs required: Sinogram, Cameras, [Attenuation], [PointSpreadFunction], [EnableGpu], [Background], [BackgroundAttenuation], [TruncateNegativeValues]");
   }

   mxClassID cid_sino  = mxGetClassID(prhs[0]);
   int       dim_sino  = mxGetNumberOfDimensions(prhs[0]); 
   
   mxClassID cid_cameras = mxGetClassID(prhs[1]);
   int       dim_cameras = mxGetNumberOfDimensions(prhs[1]);

   mxClassID cid_attenuation = mxGetClassID(prhs[2]);
   int       dim_attenuation = mxGetNumberOfDimensions(prhs[2]);

   mxClassID cid_psf = mxGetClassID(prhs[3]);
   int       dim_psf = mxGetNumberOfDimensions(prhs[3]);   

   int sino_size[3];     // Size of input sinogram matrix
   int psf_size[3];      // Size of input psf matrix
   int bkpr_size[3];     // Size of output backprojection matrix
   int cameras_size[2];  // Size of cameras matrix (can be (nx3) or (1xn))
   int enable_gpu = 0;      // Flag that enables(1)/disables(0) GPU acceleration
   float background = 0; // Background (for rotation in the backprojection)
   float background_attenuation=0; // Attenuation background (when rotating attenuation)
   int no_psf = 0;         // This flag goes high if psf input parameter from the Matlab function is a scalar -> no psf
   int no_attenuation = 0; // This flag goes high if attenuation is a scalar -> no attenuation 
   int truncate_negative_values = 1; // Set this flag to 1 in order to truncate the results to 0 if negative. Set it to 0 to disable truncation. 

   int status = 1;       // Return status: 0 if succesful

   /* The inputs must be noncomplex single floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Sinogram' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Cameras' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Attenuation' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Psf' must be noncomplex single or double.");

   /* Check if size of cameras matrix is correct */     
   cameras_size[0] = mxGetDimensions(prhs[1])[0];
   cameras_size[1] = mxGetDimensions(prhs[1])[1];
   if (!(cameras_size[1] == 3 || cameras_size[1] == 1 || cameras_size[1] == 6))
      mexErrMsgTxt("Cameras must be of size [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6]");     
   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: Backprojection");
   } 

   /* Check is attenuation is a scalar, in that case do not apply any attenuation */
   if (nrhs<3)  //no attenuation parameter
       no_attenuation = 1;
   else
       {
       cid_attenuation = mxGetClassID(prhs[2]);
       dim_attenuation = mxGetNumberOfDimensions(prhs[2]);  
       int all_one=1;
       for (int i=0; i<dim_attenuation; i++)
           {
           if (mxGetDimensions(prhs[2])[i]!=1)
              all_one = 0;
           }
       if (all_one)
           no_attenuation = 1;  //attenuation parameter is a scalar
       }

   /* Check is psf is a scalar, in that case do not apply any psf */
   if (nrhs<4)  //no psf parameter
       no_psf = 1;
   else
       {
       cid_psf = mxGetClassID(prhs[3]);
       dim_psf = mxGetNumberOfDimensions(prhs[3]);  
       int all_one=1;
       for (int i=0; i<dim_psf; i++)
           if (mxGetDimensions(prhs[3])[i]!=1)
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
            mexErrMsgTxt("For 2D reconstructions use sinograms of size [Nx1xN_projections]. E.g. sinogram = reshape(sinogram,N,1,N_projections); ");
      /* Check consistency of input if 3D (and create size of sinogram for return */
      case 3:
           if (!no_psf)
               if (dim_psf != 3)
                   mexErrMsgTxt("Dimension of Point Spread Function matrix must match the dimension of Sinogram (ddpsf).");           
           sino_size[0] = mxGetDimensions(prhs[0])[0];
           sino_size[1] = mxGetDimensions(prhs[0])[1];
           sino_size[2] = mxGetDimensions(prhs[0])[2];
//           if (sino_size[0]<2 || sino_size[1]<2)
//               mexErrMsgTxt("Size of Activity matrix must be greater then 2");           
           if (sino_size[2] != cameras_size[0])
                mexErrMsgTxt("Number of cameras must match Sinogram size");           
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[3])[0];
               psf_size[1] = mxGetDimensions(prhs[3])[1];
               psf_size[2] = mxGetDimensions(prhs[3])[2];
               if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=sino_size[0])
                    mexErrMsgTxt("Point Spread Function must be of size hxkxN for Activity of size NxmxN; h,k odd.");
               }
           bkpr_size[0] = sino_size[0];
           bkpr_size[1] = sino_size[1];
           bkpr_size[2] = sino_size[0];
           if(!no_attenuation)
               if(bkpr_size[0] != mxGetDimensions(prhs[2])[0] || bkpr_size[1] != mxGetDimensions(prhs[2])[1] || bkpr_size[2] != mxGetDimensions(prhs[2])[2])
                   mexErrMsgTxt("Attenuation must be of the same size of Activity");
           break;        
      default:
           mexErrMsgTxt("Activity must be either 2D or 3D.");
           break;
      }

   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>=5)
       enable_gpu = (int) (mxGetScalar(prhs[4]));

   /* Check if activity is multiple of ET_BLOCK_SIZE */
 //  if (enable_gpu) {
 //      if (!et_array_is_block_multiple(bkpr_size[0]) || !et_array_is_block_multiple(bkpr_size[1])) {
 //          char msg[100];
 //          sprintf(msg,"With GPU enabled, size of activity must be a multiple of %d",et_array_get_block_size());
 //          mexErrMsgTxt(msg);
 //          }
 //      }

   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=6)
       background = (mxGetScalar(prhs[5]));

   /* Check if BackgroundAttenuation is specified */
   background_attenuation = 0.0f;
   if (nrhs>=7)
       background_attenuation = (mxGetScalar(prhs[6]));

   /* Check if 'truncate negative values' flag has been specified */
   if (nrhs>=8)
       truncate_negative_values = (mxGetScalar(prhs[7])); 

   /* Extract pointers to input matrices */
   float *sinogram_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       sinogram_ptr = (float *) (mxGetData(prhs[0]));
   else
   {
       double *sinogram_ptr_d = (double *) (mxGetData(prhs[0]));
       sinogram_ptr = (float*) malloc(sino_size[0]*sino_size[1]*sino_size[2]*sizeof(float));
       for (int i=0; i<sino_size[0]*sino_size[1]*sino_size[2];i++)
           sinogram_ptr[i] = sinogram_ptr_d[i];
   }

   float *cameras_ptr;
   if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
       cameras_ptr = (float *) (mxGetData(prhs[1]));
   else
   {
       double *cameras_ptr_d = (double *) (mxGetData(prhs[1]));
       cameras_ptr = (float*) malloc(cameras_size[0]*cameras_size[1]*sizeof(float));
       for (int i=0; i<cameras_size[0]*cameras_size[1];i++)
           cameras_ptr[i] = cameras_ptr_d[i];  
   }

   int attenuation_size[3];
   float *attenuation_ptr=NULL;
   if(no_attenuation)
       {
       attenuation_size[0]=0; attenuation_size[1]=0; attenuation_size[2]=0;
       }
   else
       {
       attenuation_size[0] = bkpr_size[0];
       attenuation_size[1] = bkpr_size[1];
       attenuation_size[2] = bkpr_size[2];
       if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
           attenuation_ptr = (float *) (mxGetData(prhs[2]));
       else
           {
           double *attenuation_ptr_d = (double *) (mxGetData(prhs[2]));
           attenuation_ptr = (float *) malloc(attenuation_size[0]*attenuation_size[1]*attenuation_size[2]*sizeof(float));
           for (int i=0; i<attenuation_size[0]*attenuation_size[1]*attenuation_size[2];i++)
               attenuation_ptr[i] = attenuation_ptr_d[i];  
           }
       }

   float *psf_ptr=NULL;
   if (!no_psf)
       {
       if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
           psf_ptr = (float *) (mxGetData(prhs[3]));
       else
           {
           double *psf_ptr_d = (double *) (mxGetData(prhs[3]));
           psf_ptr = (float *) malloc(psf_size[0]*psf_size[1]*psf_size[2]*sizeof(float));
           for (int i=0; i<psf_size[0]*psf_size[1]*psf_size[2];i++)
               psf_ptr[i] = psf_ptr_d[i];  
           }
       }
           
   /* Allocate backprojection matrix */   
   int dim_bkpr;   
   mwSize mw_bkpr_size[3];
   
   dim_bkpr = dim_sino;
   mw_bkpr_size[0] = (mwSize)bkpr_size[0];
   mw_bkpr_size[1] = (mwSize)bkpr_size[1]; 
   mw_bkpr_size[2] = (mwSize)bkpr_size[2]; 

   plhs[0] =  mxCreateNumericArray(dim_bkpr, mw_bkpr_size, mxSINGLE_CLASS, mxREAL);
   float *bkpr_ptr = (float *)(mxGetData(plhs[0]));

   /* Perform backprojection */
   status = SPECT_backproject_parallelholes(sinogram_ptr, sino_size, bkpr_ptr, bkpr_size, cameras_ptr, cameras_size, psf_ptr, psf_size, attenuation_ptr, attenuation_size, &background, &background_attenuation, &enable_gpu, &truncate_negative_values);

   /* Shutdown */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(sinogram_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(cameras_ptr);
   if((!no_attenuation) && (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)) free(attenuation_ptr);
   if ((!no_psf) && (mxGetClassID(prhs[3]) != mxSINGLE_CLASS)) free(psf_ptr);
   
   /* Return */
   if (status != niftyrec_success)
   	mexErrMsgTxt(niftyrec_error_msg(status));
   return;
}


