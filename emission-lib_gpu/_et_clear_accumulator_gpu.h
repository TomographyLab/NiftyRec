/*
 *  _et_clear_accumulator_gpu.h
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

void et_clear_accumulator_gpu(float **d_accumulator, nifti_image *accumulator);
