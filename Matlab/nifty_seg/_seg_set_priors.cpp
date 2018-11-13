/*
 *  _seg_set_priors.cpp
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

#define MAX_CLASSES 30

/*###############################################################################################*/
/* Matlab extensions */
/*###############################################################################################*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int status = 0;
   int n_classes = nrhs;
   char *filenames[MAX_CLASSES];
   int buflen;

   /* Check for proper number of arguments. */
   if (n_classes < 2)
      mexErrMsgTxt("input: one file name per prior:  set_set_priors(filename1, filename2, filename3, ..)");

   /*  */
   for (int file_n = 0; file_n<n_classes; file_n++)
   {
      /* Input must be a string. */
      if (mxIsChar(prhs[file_n]) != 1)
          mexErrMsgTxt("input: one file name per prior:  set_set_priors('filename1.nii', 'filename2.nii', 'filename3.nii', ..)");

      /* Get the length of the input string. */
      buflen = (mxGetM(prhs[file_n]) * mxGetN(prhs[file_n])) + 1;

      /* Allocate memory for input and output strings. */ 
      filenames[file_n]  = (char*) calloc(buflen, sizeof(char));

      /* Copy the string data from prhs[0] into the C buffer. */
      status = mxGetString(prhs[file_n], filenames[file_n], buflen);
      if (status != 0) 
          mexWarnMsgTxt("Not enough space. String is truncated.");
   }

   status = seg_array_set_priors_fromfiles(n_classes, filenames);

   for (int file_n = 0; file_n<n_classes; file_n++)
   {
      free(filenames[file_n]);
   }

   /* Return */
   if (status != 0)
   	mexErrMsgTxt("Error while loading NiftySeg priors.");
   return;
}


