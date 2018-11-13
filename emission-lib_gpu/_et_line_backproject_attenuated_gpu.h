/*
 *  _et_line_backproject_attenuated_gpu.h 
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
#include <_reg_blocksize_gpu.h>
#include "nifti1_io.h"

#define BLOCK 256

void et_line_backproject_attenuated_gpu(float **d_sinogram, float **d_backprojection, float **d_attenuation, int cam, nifti_image *backprojection);
