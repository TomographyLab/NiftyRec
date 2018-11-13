/*
 *  _et_variance_centre_approx_mex.cpp
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
   if (!(nrhs==2 || nrhs==3 || nrhs==4 || nrhs==5 || nrhs==6 || nrhs==7)){
      mexErrMsgTxt("2, 3, 4, 5, 6 or 7 inputs required: Activity, Cameras, Attenuation, PointSpreadFunction, EnableGPU, Background, BackgroundAttenuation");
   }

   mxClassID cid_activity  = mxGetClassID(prhs[0]);
   int       dim_activity  = mxGetNumberOfDimensions(prhs[0]); 

   mxClassID cid_cameras   = mxGetClassID(prhs[1]);
   int       dim_cameras   = mxGetNumberOfDimensions(prhs[1]); 

   mxClassID cid_psf;
   int       dim_psf          = 0;

   mxClassID cid_attenuation;
   int       dim_attenuation  = 0;

   int activity_size[3];  //size of input activity matrix.
   int attenuation_size[3];//size of attenuation matrix.
   int psf_size[3];       //size of input psf matrix.
   int sino_size[3];      //size of output sinogram matrix.
   int cameras_size[2];   //size of input cameras matrix.
   float background;      // Background value, defaults to 0. (It is used to fill empty areas when images are rotated)
   float background_attenuation; // Attenuation background value, defaults to 0.
   int enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int no_psf = 0;        // This flag goes high if psf input parameter from the Matlab function is a scalar -> no psf
   int no_attenuation = 0;// This flag goes high if an attenuation image is given and it is not a scalar. 

   int N;                 // Size of activity: [N,N,m].
   int m;                 // Size of activity: [N,N,m].
   int status = 1;        // Return 0 if everython ok, 1 if errors.

   int temp;

   /* The inputs must be noncomplex single floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Activity' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Cameras' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Attenuation' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Psf' must be noncomplex single or double.");

   /* Check if size of cameras matrix is correct */
   cameras_size[0] = mxGetDimensions(prhs[1])[0];
   cameras_size[1] = mxGetDimensions(prhs[1])[1];
   if (!(cameras_size[1] == 3 || cameras_size[1] == 1))
      mexErrMsgTxt("Cameras must be of size [n_cameras x 1] or [n_cameras x 3]");        
   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: Sinogram");
   }     
   /* Check is psf is a scalar, in that case do not apply any psf */
   if (nrhs<4)  //no psf parameter
       no_psf = 1;
   else
       {
       cid_psf = mxGetClassID(prhs[3]);
       dim_psf = mxGetNumberOfDimensions(prhs[3]);  
       temp=1;
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
       cid_attenuation = mxGetClassID(prhs[2]);
       dim_attenuation = mxGetNumberOfDimensions(prhs[2]);  
       temp=1;
       for (int i=0; i<dim_attenuation; i++)
           {
           if (mxGetDimensions(prhs[2])[i]!=1)
              temp = 0;
           }
       if (temp)
           no_attenuation = 1;  //attenuation parameter is a scalar
       }

   /* Check consistency of input (and create size of sinogram for return */
   switch(dim_activity)
      {
      /* Check consistency of input if 2D (and create size of sinogram for return) */
      case 2:
           if (dim_psf != 2 && !no_psf)
               mexErrMsgTxt("Size of Point Spread Function matrix must match the size of Activity (ddpsf).");        
           activity_size[0] = mxGetDimensions(prhs[0])[0];
           activity_size[1] = mxGetDimensions(prhs[0])[1];
           activity_size[2] = 1;
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[3])[0];
               psf_size[1] = mxGetDimensions(prhs[3])[1];
               psf_size[2] = 1;
               if ((psf_size[0]%2!=1) || (psf_size[1]!=activity_size[0]))
                   mexErrMsgTxt("For activity of size (NxN) Point Spread Function must be of size hxN; h odd");
               }
           N = activity_size[0];
           if (N<2)
               mexErrMsgTxt("Size of Activity matrix must be greater then 2");
           if (activity_size[1]!=N)
               mexErrMsgTxt("2D Activity must be square (3D activity can be of size NxNxm)");
           sino_size[0] = N;
           sino_size[1] = cameras_size[0];
           sino_size[2] = 1;
           break;
           if(!no_attenuation)
               if(mxGetDimensions(prhs[0])[0] != mxGetDimensions(prhs[2])[0] || mxGetDimensions(prhs[0])[1] != mxGetDimensions(prhs[2])[1])
                   mexErrMsgTxt("Attenuation must be of the same size of Activity");
      /* Check consistency of input if 3D (and create size of sinogram for return */
      case 3:
           if (!no_psf)
               if (dim_psf != 3 && !no_psf)
                   mexErrMsgTxt("Size of Point Spread Function matrix must match the size of Activity (ddpsf).");           
           activity_size[0] = mxGetDimensions(prhs[0])[0];
           activity_size[1] = mxGetDimensions(prhs[0])[1];
           activity_size[2] = mxGetDimensions(prhs[0])[2];
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[3])[0];
               psf_size[1] = mxGetDimensions(prhs[3])[1];
               psf_size[2] = mxGetDimensions(prhs[3])[2];
               if ((psf_size[0]%2!=1) || ((psf_size[1]%2!=1)) || psf_size[2]!=activity_size[2])
                   mexErrMsgTxt("Point Spread Function must be of size hxkxm for Activity of size NxNxm; h,k odd");
               }
           N = activity_size[0];
           m = activity_size[2];
           if (activity_size[0]<2 || activity_size[1]<2 || activity_size[2]<2)
               mexErrMsgTxt("Size of Activity matrix must be greater then 2");           
           if (activity_size[1]!=activity_size[0])
               mexErrMsgTxt("2D Activity must be square (3D activity can be of size NxNxm)");
           sino_size[0] = activity_size[0];
           sino_size[1] = activity_size[2];
           sino_size[2] = cameras_size[0];
           if(!no_attenuation)
               if(mxGetDimensions(prhs[0])[0] != mxGetDimensions(prhs[2])[0] || mxGetDimensions(prhs[0])[1] != mxGetDimensions(prhs[2])[1] || mxGetDimensions(prhs[0])[2] != mxGetDimensions(prhs[2])[2])
                   mexErrMsgTxt("Attenuation must be of the same size of Activity");
           break;        
      default:
           mexErrMsgTxt("Activity must be either 2D or 3D.");
           break;
      }

   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>=4)
       enable_gpu = (int) (mxGetScalar(prhs[4]));
  
   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>=5)
       background = (mxGetScalar(prhs[5])); 

   /* Check if attenuation background is specified */
   background_attenuation = 0.0f;
   if (nrhs>=6)
       background_attenuation = (mxGetScalar(prhs[6])); 

   /* Extract pointers to input matrices (and eventually convert double to float)*/
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

   if(no_attenuation)
       {
       attenuation_size[0]=0; attenuation_size[1]=0; attenuation_size[2]=0;
       }
   else
       {
       attenuation_size[0]=activity_size[0]; attenuation_size[1]=activity_size[1]; attenuation_size[2]=activity_size[2];
       }

   float *attenuation_ptr;
   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       attenuation_ptr = (float *) (mxGetData(prhs[2]));
   else
   {
       double *attenuation_ptr_d = (double *) (mxGetData(prhs[2]));
       attenuation_ptr = (float*) malloc(attenuation_size[0]*attenuation_size[1]*attenuation_size[2]*sizeof(float));
       for (int i=0; i<attenuation_size[0]*attenuation_size[1]*attenuation_size[2];i++)
           attenuation_ptr[i] = attenuation_ptr_d[i];
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

   float *psf_ptr;
   if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
       psf_ptr = (float *) (mxGetData(prhs[3]));
   else
   {
       double *psf_ptr_d = (double *) (mxGetData(prhs[3]));
       psf_ptr = (float *) malloc(psf_size[0]*psf_size[1]*psf_size[2]*sizeof(float));
       for (int i=0; i<psf_size[0]*psf_size[1]*psf_size[2];i++)
           psf_ptr[i] = psf_ptr_d[i];  
   }

   /* Check if rotations are not only along Z axis: in this case activity must be a cube */
//   if (dim_activity != 3 || activity_size[2] != activity_size[0])
//       if (cameras_size[1]==3)
//           for (int cam=0; cam<cameras_size[0]; cam++)
//               if (fabs(cameras_ptr[1*cameras_size[0]+cam])>eps || fabs(cameras_ptr[2*cameras_size[0]+cam])>eps)
//                   mexErrMsgTxt("At least one of the cameras has multiple axis of rotation, in this case Activity must be a cube (N,N,N)");

   /* Allocate sinogram matrix */   
   int dim_sino;   
   mwSize mw_sino_size[3];
   
   dim_sino = dim_activity;
   mw_sino_size[0] = (mwSize)sino_size[0];
   mw_sino_size[1] = (mwSize)sino_size[1]; 
   mw_sino_size[2] = (mwSize)sino_size[2]; 

   plhs[0] =  mxCreateNumericArray(dim_sino, mw_sino_size, mxSINGLE_CLASS, mxREAL);
   float *sinogram_ptr = (float *)(mxGetData(plhs[0]));

   /* Perform projection */
//   status = et_array_approx_variance_centre(activity_ptr, activity_size, sinogram_ptr, sino_size, cameras_ptr, cameras_size, psf_ptr, psf_size, attenuation_ptr, attenuation_size, background, background_attenuation, enable_gpu);

   /* Shutdown */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(activity_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(cameras_ptr);
   if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) free(attenuation_ptr);
   if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) free(psf_ptr);

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while performing projection.");
   return;
}


