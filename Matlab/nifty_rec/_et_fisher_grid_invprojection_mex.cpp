/*
 *  _et_fisher_grid_invprojection_mex.cpp
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
   if (!(nrhs==3 || nrhs==4 || nrhs==5 || nrhs==6 || nrhs==7 || nrhs==8 || nrhs==9)){
      mexErrMsgTxt("3, 4, 5, 6, 7 or 8 inputs required: Inverse (element by element) of the Projection, Cameras, Grid, Attenuation, PointSpreadFunction, EnableGPU, Epsilon, Background, BackgroundAttenuation");
   }

   mxClassID cid_projection  = mxGetClassID(prhs[0]);
   int       dim_projection  = mxGetNumberOfDimensions(prhs[0]); 

   mxClassID cid_cameras   = mxGetClassID(prhs[1]);
   int       dim_cameras   = mxGetNumberOfDimensions(prhs[1]); 

   mxClassID cid_psf;
   int       dim_psf          = 0;

   mxClassID cid_attenuation;
   int       dim_attenuation  = 0;

   int projection_size[3];  //size of input projection matrix.
   int attenuation_size[3];//size of attenuation matrix.
   int psf_size[3];       //size of input psf matrix.
   int bkpr_size[3];     // Size of output backprojection matrix
   int cameras_size[2];   //size of input cameras matrix.
   float background;      // Background value, defaults to 0. (It is used to fill empty areas when images are rotated)
   float background_attenuation; // Attenuation background value, defaults to 0.
   float epsilon;
   int enable_gpu;        // Flag for GPU Acceleration: 1 to enable, 0 to disable.
   int no_psf = 0;        // This flag goes high if psf input parameter from the Matlab function is a scalar -> no psf
   int no_attenuation = 0;// This flag goes high if an attenuation image is given and it is not a scalar. 

//   int N;                 // Size of projection
//   int m;                 // Size of projection
   int status = 1;        // Return 0 if everython ok, 1 if errors.

   int temp;

   /* The inputs must be noncomplex single floating point matrices */
   if ( (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Projection' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Cameras' must be noncomplex single or double.");
   if ( (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Grid' must be noncomplex single or double.");
   if (nrhs>=4)
       if ( (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Attenuation' must be noncomplex single or double.");
   if (nrhs>=5)
       if ( (mxGetClassID(prhs[4]) != mxSINGLE_CLASS) && (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) ) mexErrMsgTxt("'Psf' must be noncomplex single or double.");

   /* Check if size of cameras matrix is correct */
   cameras_size[0] = mxGetDimensions(prhs[1])[0];
   cameras_size[1] = mxGetDimensions(prhs[1])[1];
   if (!(cameras_size[1] == 3 || cameras_size[1] == 1))
      mexErrMsgTxt("Cameras must be of size [n_cameras x 1] or [n_cameras x 3]");        
   /* Check if number of output arguments is correct */
   if (!(nlhs == 1 || nlhs == 2)){
      mexErrMsgTxt("One or two outputs: [Fisher Information Matrix, Quadratic Prior Fisher Information Matrix]");
   }     

   /* Check is psf is a scalar, in that case do not apply any psf */
   if (nrhs<5)  //no psf parameter
       no_psf = 1;
   else
       {
       cid_psf = mxGetClassID(prhs[4]);
       dim_psf = mxGetNumberOfDimensions(prhs[4]);  
       temp=1;
       for (int i=0; i<dim_psf; i++)
           if (mxGetDimensions(prhs[4])[i]!=1)
              temp = 0;
       if (temp)
           no_psf = 1;  //psf parameter is a scalar
       }

   /* Require PSF */
   if (no_psf)
       mexErrMsgTxt("PSF must be specified in order to compute the Fisher Information Matrix.");

   /* Check is attenuation is a scalar, in that case do not apply any attenuation */
   if (nrhs<4)  //no attenuation parameter
       no_attenuation = 1;
   else
       {
       cid_attenuation = mxGetClassID(prhs[3]);
       dim_attenuation = mxGetNumberOfDimensions(prhs[3]);  
       temp=1;
       for (int i=0; i<dim_attenuation; i++)
           {
           if (mxGetDimensions(prhs[3])[i]!=1)
              temp = 0;
           }
       if (temp)
           no_attenuation = 1;  //attenuation parameter is a scalar
       }

   /* Check consistency of input (and create size of sinogram for return */
   switch(dim_projection)
      {
      /* Check consistency of input if 2D (and create size of sinogram for return) */
      case 2:
           if (!no_psf)
               if (dim_psf != 2)
                   mexErrMsgTxt("Dimension of Point Spread Function matrix must match the dimension of Projection (ddpsf).");        
           projection_size[0] = mxGetDimensions(prhs[0])[0];
           projection_size[1] = mxGetDimensions(prhs[0])[1];
           projection_size[2] = 1;
//           if (projection_size[0]<2)
//               mexErrMsgTxt("Size of Activity matrix must be greater then 2");
           if (projection_size[1] != cameras_size[0])
               mexErrMsgTxt("Number of cameras must match Sinogram size");    
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[4])[0];
               psf_size[1] = mxGetDimensions(prhs[4])[1];
               psf_size[2] = 1;
               if (psf_size[0]%2!=1 || psf_size[1]!=projection_size[0])
                   mexErrMsgTxt("Point Spread Function must be of size hxN for Backprojection of size Nxn; h odd.");
               }
           bkpr_size[0] = projection_size[0];
           bkpr_size[1] = projection_size[0];
           bkpr_size[2] = 1;
           if(!no_attenuation)
               if(bkpr_size[0] != mxGetDimensions(prhs[3])[0] || bkpr_size[1] != mxGetDimensions(prhs[3])[1])
                   mexErrMsgTxt("Attenuation must be of the same size of Activity");
           if(bkpr_size[0] != mxGetDimensions(prhs[2])[0] || bkpr_size[1] != mxGetDimensions(prhs[2])[1])
                   mexErrMsgTxt("Grid must be of the same size of Activity");
           break;
      /* Check consistency of input if 3D (and create size of sinogram for return */
      case 3:
           if (!no_psf)
               if (dim_psf != 3)
                   mexErrMsgTxt("Dimension of Point Spread Function matrix must match the dimension of Projection (ddpsf).");           
           projection_size[0] = mxGetDimensions(prhs[0])[0];
           projection_size[1] = mxGetDimensions(prhs[0])[1];
           projection_size[2] = mxGetDimensions(prhs[0])[2];
           if (projection_size[0]<2 || projection_size[1]<2)
               mexErrMsgTxt("Size of Activity matrix must be greater then 2");           
           if (projection_size[2] != cameras_size[0])
                mexErrMsgTxt("Number of cameras must match projection size");           
           if (!no_psf)
               {
               psf_size[0] = mxGetDimensions(prhs[4])[0];
               psf_size[1] = mxGetDimensions(prhs[4])[1];
               psf_size[2] = mxGetDimensions(prhs[4])[2];
               if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=projection_size[0])
                    mexErrMsgTxt("Point Spread Function must be of size hxkxN for Activity of size NxmxN; h,k odd.");
               }
           bkpr_size[0] = projection_size[0];
           bkpr_size[1] = projection_size[0];
           bkpr_size[2] = projection_size[1];
           if(!no_attenuation)
               if(bkpr_size[0] != mxGetDimensions(prhs[3])[0] || bkpr_size[1] != mxGetDimensions(prhs[3])[1] || bkpr_size[2] != mxGetDimensions(prhs[3])[2])
                   mexErrMsgTxt("Attenuation must be of the same size of Activity");
           if(bkpr_size[0] != mxGetDimensions(prhs[2])[0] || bkpr_size[1] != mxGetDimensions(prhs[2])[1] || bkpr_size[2] != mxGetDimensions(prhs[2])[2])
                   mexErrMsgTxt("Grid must be of the same size of Activity");

           break;        
      default:
           mexErrMsgTxt("Activity must be either 2D or 3D.");
           break;
      }

   /* Check if activity is multiple of ET_BLOCK_SIZE */
//   if (!et_is_block_multiple(bkpr_size[0]) || !et_is_block_multiple(bkpr_size[1]))
//       mexErrMsgTxt("Size of activity must be a multiple of 64");

   /* Check if EnableGPU is specified */
   enable_gpu = 0;
   if (nrhs>5)
       enable_gpu = (int) (mxGetScalar(prhs[5]));

   /* Check if Epsilon is specified */
   epsilon = 0.0f;
   if (nrhs>6)
       epsilon = (mxGetScalar(prhs[6])); 
  
   /* Check if Background is specified */
   background = 0.0f;
   if (nrhs>7)
       background = (mxGetScalar(prhs[7])); 

   /* Check if attenuation background is specified */
   background_attenuation = 0.0f;
   if (nrhs>8)
       background_attenuation = (mxGetScalar(prhs[8])); 

   /* Extract pointers to input matrices (and eventually convert double to float)*/
   float *projection_ptr;
   if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS)
       projection_ptr = (float *) (mxGetData(prhs[0]));
   else
   {
       double *projection_ptr_d = (double *) (mxGetData(prhs[0]));
       projection_ptr = (float*) malloc(projection_size[0]*projection_size[1]*projection_size[2]*sizeof(float));
       for (int i=0; i<projection_size[0]*projection_size[1]*projection_size[2];i++)
           projection_ptr[i] = projection_ptr_d[i];
   }

   float *attenuation_ptr=NULL;
   if(no_attenuation)
       {
       attenuation_size[0]=0; attenuation_size[1]=0; attenuation_size[2]=0;
       }
   else
       {
       attenuation_size[0]=bkpr_size[0]; attenuation_size[1]=bkpr_size[1]; attenuation_size[2]=bkpr_size[2];

       if (mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
           attenuation_ptr = (float *) (mxGetData(prhs[3]));
       else
           {
           double *attenuation_ptr_d = (double *) (mxGetData(prhs[3]));
           attenuation_ptr = (float*) malloc(attenuation_size[0]*attenuation_size[1]*attenuation_size[2]*sizeof(float));
           for (int i=0; i<attenuation_size[0]*attenuation_size[1]*attenuation_size[2];i++)
               attenuation_ptr[i] = attenuation_ptr_d[i];
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

   float *psf_ptr;
   if (no_psf) 
       {
       psf_size[0]=0; psf_size[1]=0; psf_size[2]=0;  
       }
   else
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

   float *grid_ptr;
   if (mxGetClassID(prhs[2]) == mxSINGLE_CLASS)
       grid_ptr = (float *) (mxGetData(prhs[2]));
   else
   {
       double *grid_ptr_d = (double *) (mxGetData(prhs[2]));
       grid_ptr = (float*) malloc(bkpr_size[0]*bkpr_size[1]*bkpr_size[2]*sizeof(float));
       for (int i=0; i<bkpr_size[0]*bkpr_size[1]*bkpr_size[2];i++)
           grid_ptr[i] = grid_ptr_d[i];
   }

   /* Allocate Fisher information matrix */   
   mwSize mw_fisher_size[2];
   int grid_elements=0;
   for (int i=0; i<bkpr_size[0]*bkpr_size[1]*bkpr_size[2];i++)
       if (grid_ptr[i]!=0)
           grid_elements++;
   mw_fisher_size[0] = (mwSize)grid_elements;
   mw_fisher_size[1] = (mwSize)grid_elements; 
   int fisher_size[2];
   fisher_size[0] = grid_elements;
   fisher_size[1] = grid_elements;

   fprintf(stderr,"Size Fisher: %d \n\n",grid_elements);
   plhs[0] =  mxCreateNumericArray(2, mw_fisher_size, mxSINGLE_CLASS, mxREAL);
   float *fisher_ptr = (float *)(mxGetData(plhs[0]));

   /* Allocate Fisher Information Matrix of the prior */
   float *fisher_prior_ptr = NULL;
   if (nlhs>=2)
       {
       plhs[1] =  mxCreateNumericArray(2, mw_fisher_size, mxSINGLE_CLASS, mxREAL);
       fisher_prior_ptr = (float *)(mxGetData(plhs[1]));
       }

   /* Calculate Fisher Information matrix */
   if (no_psf)
       mexErrMsgTxt("PSF must be specified in order to compute the Fisher Information Matrix.");
   status = et_array_fisher_grid_projection(projection_ptr, projection_size, bkpr_size, cameras_ptr, cameras_size, psf_ptr, psf_size, grid_ptr, fisher_ptr, fisher_prior_ptr, fisher_size, attenuation_ptr, attenuation_size, epsilon, background, background_attenuation, enable_gpu);

   /* Shutdown */
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS) free(projection_ptr);
   if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) free(cameras_ptr);
   if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS) free(grid_ptr);
   if(!no_attenuation) if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) free(attenuation_ptr);
   if (!no_psf) if (mxGetClassID(prhs[4]) != mxSINGLE_CLASS) free(psf_ptr);

   /* Return */
   if (status == ET_ERROR_BADGRID)
        mexErrMsgTxt("Error while calculating the Fisher Information matrix from projection: the N non-zero elements of the grid must be unique integers in the range from 1 to N-1.");
   if (status != 0)
   	mexErrMsgTxt("Error while calculating the Fisher Information matrix from projection.");
   return;
}



