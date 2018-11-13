/*
 *  _et_line_backproject_attenuated.h
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

void et_line_backproject_attenuated(nifti_image *sinogramImage, nifti_image *temp_backprojectionImage, nifti_image *attenuationImage, int cam);

