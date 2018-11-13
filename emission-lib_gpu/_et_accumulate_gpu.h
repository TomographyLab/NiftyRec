/*
 *  _et_accumulate_gpu.h
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

#define BLOCK 256

/*  d_B = d_B + d_A  */
void et_accumulate_gpu(float **d_A, float **d_B, nifti_image *n_image);
