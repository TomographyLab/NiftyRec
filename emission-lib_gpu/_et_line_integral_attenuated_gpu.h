/*
 *  _et_line_integrals_attenuated_gpu.h
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
#include "_reg_blocksize_gpu.h"
#include "nifti1_io.h"

void et_line_integral_attenuated_gpu(float *d_activity, float *d_attenuation, float *d_sinogram, float *d_background_image, float *d_partialsum, int cam, nifti_image *, float background_activity);

