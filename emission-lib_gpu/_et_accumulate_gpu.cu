/*
 *  _et_accumulate_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_accumulate_gpu.h"
#include "_et_accumulate_gpu_kernels.cu"

/*  d_B = d_B + d_A  */
void et_accumulate_gpu(float **d_A, float **d_B, nifti_image *image)
{
	int3 image_size = make_int3(image->nx,image->ny,image->nz);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_image_size,&image_size,sizeof(int3)));

	const unsigned int grid = (unsigned int)ceil(image->nx*image->ny*image->nz/(float)BLOCK);
	dim3 B(BLOCK,1,1);
	dim3 G(grid,1,1);

	et_accumulate_gpu_kernel <<<G,B>>> (*d_A,*d_B);

	CUDA_SAFE_CALL(cudaThreadSynchronize());
}			

