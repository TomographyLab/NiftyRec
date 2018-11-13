/*
 *  _et_reorder_activity_out.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing
 *  UCL - University College London 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <mex.h>

/*###############################################################################################*/
/* Matlab extension */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Check for proper number of arguments. */
   if (nrhs!=1){
      mexErrMsgTxt("One input: 2D or 3D Activity matrix");
   }

   mxClassID cid_activity  = mxGetClassID(prhs[0]);
   int dim_activity  = mxGetNumberOfDimensions(prhs[0]); 
   int activity_size[3];  //size of input activity matrix.

   /* Extract pointer to input matrix */
   float *activity_ptr = (float *) (mxGetData(prhs[0]));
   if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
       mexErrMsgTxt("Input must be single");

   /* Check if number of output arguments is correct */
   if (nlhs != 1){
      mexErrMsgTxt("One output: 2D or 3D Activity matrix");
   }    

   /* Allocate output */
   switch(dim_activity)
   {
      case 2:
           activity_size[0] = mxGetDimensions(prhs[0])[0];
           activity_size[1] = mxGetDimensions(prhs[0])[1];
           activity_size[2] = 1;
           break;
      case 3:
           activity_size[0] = mxGetDimensions(prhs[0])[0];
           activity_size[1] = mxGetDimensions(prhs[0])[1];
           activity_size[2] = mxGetDimensions(prhs[0])[2];;
           break;    
   }
 
   mwSize mw_activity_size[3];
   mw_activity_size[0] = (mwSize)activity_size[1];  //x'=y
   mw_activity_size[1] = (mwSize)activity_size[2];  //y'=z
   mw_activity_size[2] = (mwSize)activity_size[0];  //z'=x

   plhs[0] =  mxCreateNumericArray(dim_activity, mw_activity_size, mxDOUBLE_CLASS, mxREAL);
   double *activity_double = (double *)(mxGetData(plhs[0]));

   int k = 0;
   for (int j=0; j<activity_size[0]*activity_size[1]; j++)
   {
       for (int i=0; i<activity_size[2]; i++)
       {
           activity_double[j+i*activity_size[0]*activity_size[1]] = activity_ptr[k];
           k++;
       }
   }

   return;
}


