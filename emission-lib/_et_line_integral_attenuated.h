/*
 *  _et_line_integral_attenuated.h 
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nifti1_io.h"

void et_line_integral_attenuated(nifti_image *inputImage, nifti_image *attenuationImage, nifti_image *sinoImage, int cam, float background);
