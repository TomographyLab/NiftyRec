/*
 *  _et_line_backproject_attenuated_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_line_backproject_attenuated_gpu.h"
#include "_et_line_backproject_attenuated_gpu_kernels.cu"

void et_line_backproject_attenuated_gpu(float **d_sinogram, float **d_backprojection, float **d_attenuation, int cam, nifti_image *backprojection)
{
	int3 backprojection_size = make_int3(backprojection->nx,backprojection->ny,backprojection->nz);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_backprojection_size,&backprojection_size,sizeof(int3)));
	
	const unsigned int grid = (unsigned int)ceil(backprojection->nx*backprojection->ny/(float)BLOCK);
	dim3 B(BLOCK,1,1);
	dim3 G(grid,1,1);
	
	float *d_sinogram_ptr = (*d_sinogram) + cam * backprojection->nx * backprojection->ny;
	
	et_line_backproject_attenuated_gpu_kernel <<<G,B>>> (d_sinogram_ptr, *d_backprojection, *d_attenuation);
	
	CUDA_SAFE_CALL(cudaThreadSynchronize());
}



