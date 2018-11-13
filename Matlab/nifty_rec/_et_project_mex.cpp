/*
 *  _et_project_mex.cpp
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>


/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (!(nrhs==2 || nrhs==3 || nrhs==4 || nrhs==5 || nrhs==6 || nrhs==7 || nrhs==8)){
      mexErrMsgTxt("2, 3, 4, 5, 6, 7 or 8 inputs required: Activity, Cameras, [Attenuation], [PointSpreadFunction], [EnableGPU], [Background], [BackgroundAttenuation], [TruncateNegativeValues]");
   }

   mxClassID cid_activity  = mxGetClassID(prhs[0]);
   int       dim_activity  = mxGetNumberOfDimensions(prhs[0]); 

   mxClassID cid_cameras   = mxGetClassID(prhs[1]);
   int       dim_cameras   = mxGetNumberOfDimensions(prhs[1]); 

   mxClassID cid_psf;
   int       dim_psf          = 0;

   mxClassID cid_attenuation;
   int       dim_attenuation  = 0;

   int activity_size[3];   //size of input activity matrix.
   int attenuation_size[3];//size of attenuation matrix.
   int psf_size[3];        //size of input psf matrix.
   int sino_size[3];       //size of output sinogram matrix.
   int cameras_size[2];    //size of input cameras matrix.
   float background;       // Background value, defaults to 0. (It is used to fill empty areas when images are rotated)
   float background_attenuation; // Attenuation background value, defaults to 0.
   int enable_gpu;         // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int no_psf = 0;         // This flag goes high if psf input parameter from the Matlab function is a scalar -> no psf
   int no_attenuation = 0; // This flag goes high if attenuation is a scalar -> no attenuation 
   int no_activity = 0;    // This flag goes high if activity is a scalar -> no activity
   int truncate_negative_values = 1; // Set this flag to 1 in order to truncate the results to 0 if negative. Set it to 0 to disable truncation. 

   int N;                  // Size of activity: [N,N,m].
   int m;                  // Size of activity: [N,N,m].
   int status = 1;         // Return 0 if everython ok, 1 if errors.
   int temp = 1; 

   /* The inputs must be noncomplex single floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Activity' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Cameras' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Attenuation' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Psf' must be noncomplex single or double.");

   /* Check if size of cameras matrix is correct */
   cameras_size[0] = mxGetDimensions(prhs[1])[0];
   cameras_size[1] = mxGetDimensions(prhs[1])[1];
   if (!(cameras_size[1] == 3 || cameras_size[1] == 1 || cameras_size[1] == 6))
      mexErrMsgTxt("Cameras must be of size [n_cameras x 1], [n_cameras x 3] or [n_cameras x 6]");        
   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: Sinogram");
   }     
   /* Check if activity is a matrix or scalar - if scalar, no activity */
   temp=1;
   for (int i=0; i<dim_activity; i++)
           {
           if (mxGetDimensions(prhs[0])[i]!=1)
              temp = 0;
           }
       if (temp)
           no_activity = 1;  //activity parameter is a scalar 
   /* Check is psf is a scalar, in that case do not apply any psf */
   if (nrhs<4)  //no psf parameter
       no_psf = 1;
   else
       {
       temp=1;
       cid_psf = mxGetClassID(prhs[3]);
       dim_psf = mxGetNumberOfDimensions(prhs[3]);  
       for (int i=0; i<dim_psf; i++)
           if (mxGetDimensions(prhs[3])[i]!=1)
              temp = 0;
       if (temp)
           no_psf = 1;  //psf parameter is a scalar
       }
   if (no_psf == 1)
       {
       psf_size[0] = 0;
       psf_size[1] = 0;
       psf_size[2] = 0;
       }
   /* Check is attenuation is a scalar, in that case do not apply any attenuation */
   if (nrhs<3)  //no psf parameter
       no_attenuation = 1;
   else
       {
       temp=1;
       cid_attenuation = mxGetClassID(prhs[2]);
       dim_attenuation = mxGetNumberOfDimensions(prhs[2]);  
       for (int i=0; i<dim_attenuation; i++)
           {
           if (mxGetDimensions(prhs[2])[i]!=1)
              temp = 0;
           }
       if (temp)
           no_attenuation = 1;  //attenuation parameter is a scalar
       }

   /* Check consistency of input (and create size of sinogram for return */
   if(no_activity && no_attenuation)
       mexErrMsgTxt("At least one between 'activity' and 'attenuation' must be defined - scalar is 'not defined'. ");

   switch(dim_activity)
      {
      case 2:
            if (!no_activity)
                mexErrMsgTxt("For 2D reconstructions use images of size [Nx1xN]. E.g. activity=reshape(activity,N,1,N); ");
      /* Check consistency of input if 3D (and create size of sinogram for return */
      case 3:
           if (!no_psf)
               if (dim_psf != 3 && !no_psf)
                   mexErrMsgTxt("Size of Point Spread Function matrix must match the size of Activity (ddpsf).");  
           if (!no_activity)
           {         
               activity_size[0] = mxGetDimensions(prhs[0])[0];
               activity_size[1] = mxGetDimensions(prhs[0])[1];
               activity_size[2] = mxGetDimensions(prhs[0])[2];
           }
           else
           {         
               activity_size[0] = mxGetDimensions(prhs[2])[0];
               activity_size[1] = mxGetDimensions(prhs[2])[1];
               activity_size[2] = mxGetDimensions(prhs[2])[2];
           }
           if ((activity_size[0]!=activity_size[2]))
               mexErrMsgTxt("3D Activity must be of size [NxMxN]");
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[3])[0];
               psf_size[1] = mxGetDimensions(prhs[3])[1];
               psf_size[2] = mxGetDimensions(prhs[3])[2];
               if ((psf_size[0]%2!=1) || ((psf_size[1]%2!=1)) || psf_size[2]!=activity_size[0])
                   mexErrMsgTxt("Point Spread Function must be of size hxkxN for Activity of size NxMxN; h,k odd");
               }
           N = activity_size[0];
           m = activity_size[1];
       
           if (activity_size[2]!=activity_size[0])
               mexErrMsgTxt("2D Activity must be square (3D activity can be of size NxNxm)");
           sino_size[0] = activity_size[0];
           sino_size[1] = activity_size[1];
           sino_size[2] = cameras_size[0];
           if((!no_attenuation) && (!no_activity))
               if(mxGetDimensions(prhs[0])[0] != mxGetDimensions(prhs[2])[0] || mxGetDimensions(prhs[0])[1] != mxGetDimensions(prhs[2])[1] || mxGetDimensions(prhs[0])[2] != mxGetDimensions(prhs[2])[2])
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
//   if (enable_gpu) {
//       if (!et_array_is_block_multiple(activity_size[0]) || !et_array_is_block_multiple(activity_size[1])) {
//           char msg[100];
//           sprintf(msg,"With GPU enabled, size of activity must be a multiple of %d",et_array_get_block_size());
//           mexErrMsgTxt(msg);
//           }
//       }

   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=6)
       background = (mxGetScalar(prhs[5])); 

   /* Check if attenuation background is specified */
   background_attenuation = 0.0f;
   if (nrhs>=7)
       background_attenuation = (mxGetScalar(prhs[6])); 

   /* Check if 'truncate negative values' flag has been specified */
   if (nrhs>=8)
       truncate_negative_values = (mxGetScalar(prhs[7])); 

   /* Extract pointers to input matrices (and eventually convert double to float)*/
   float *attenuation_ptr=NULL;
   if(no_attenuation)
       {
       attenuation_size[0]=0; attenuation_size[1]=0; attenuation_size[2]=0;
       }
   else
       {
       attenuation_size[0]=activity_size[0]; attenuation_size[1]=activity_size[1]; attenuation_size[2]=activity_size[2];
       if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
           attenuation_ptr = (float *) (mxGetData(prhs[2]));
       else
           {
           double *attenuation_ptr_d = (double *) (mxGetData(prhs[2]));
           attenuation_ptr = (float*) malloc(attenuation_size[0]*attenuation_size[1]*attenuation_size[2]*sizeof(float));
           for (int i=0; i<attenuation_size[0]*attenuation_size[1]*attenuation_size[2];i++)
               attenuation_ptr[i] = attenuation_ptr_d[i];
           }
       }

   float *activity_ptr=NULL; 
   if(no_activity)
       {
       activity_size[0]=0; activity_size[1]=0; activity_size[2]=0;
       }
   else
       {
       if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
           activity_ptr = (float *) (mxGetData(prhs[0]));
       else
           {
           double *activity_ptr_d = (double *) (mxGetData(prhs[0]));
           activity_ptr = (float*) malloc(activity_size[0]*activity_size[1]*activity_size[2]*sizeof(float));
           for (int i=0; i<activity_size[0]*activity_size[1]*activity_size[2];i++)
               activity_ptr[i] = activity_ptr_d[i];
           }
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

   float *psf_ptr=NULL;
   if(!no_psf)
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

   /* Allocate sinogram matrix */   
   int dim_sino=3;   
   mwSize mw_sino_size[3];
   
   mw_sino_size[0] = (mwSize)sino_size[0];
   mw_sino_size[1] = (mwSize)sino_size[1]; 
   mw_sino_size[2] = (mwSize)sino_size[2]; 

   plhs[0] =  mxCreateNumericArray(dim_sino, mw_sino_size, mxSINGLE_CLASS, mxREAL);
   float *sinogram_ptr = (float *)(mxGetData(plhs[0]));

   /* Perform projection */
   status = SPECT_project_parallelholes(activity_ptr, activity_size, sinogram_ptr, sino_size, cameras_ptr, cameras_size, psf_ptr, psf_size, attenuation_ptr, attenuation_size, &background, &background_attenuation, &enable_gpu, &truncate_negative_values);

   /* Shutdown */
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(cameras_ptr);
   if ((!no_activity)    && (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)) free(activity_ptr);
   if ((!no_attenuation) && (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)) free(attenuation_ptr);
   if ((!no_psf)         && (mxGetClassID(prhs[3]) != mxSINGLE_CLASS)) free(psf_ptr);

   /* Return */
   if (status != niftyrec_success)
   	mexErrMsgTxt(niftyrec_error_msg(status));
   return;
}


