/*
 *  _et_clear_accumulator_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_clear_accumulator_gpu.h"

#define BLOCK 256

void et_clear_accumulator_gpu(float **d_accumulator, nifti_image *accumulator)
{
	cudaMemset((void*)*d_accumulator,0.0f,accumulator->nx*accumulator->ny*accumulator->nz*sizeof(float));
}				

