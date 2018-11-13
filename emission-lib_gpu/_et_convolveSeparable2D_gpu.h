/*
 *  _et_convolveSeparable2D_gpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <cutil_inline.h>
#include "nifti1_io.h"

#define MAX_SEPARABLE_KERNEL_RADIUS 16
#define ET_BLOCK_SIZE 8

extern "C" int convolutionRowsGPU(float *d_Dst,float *d_Src,int imageW,int imageH,int kernelRadius);
extern "C" int convolutionColumnsGPU(float *d_Dst,float *d_Src,int imageW,int imageH,int kernelRadius);

int et_convolveSeparable2D_gpu(float **d_data, int *data_size, float **d_kernel, int *kernel_size, float **d_result);


