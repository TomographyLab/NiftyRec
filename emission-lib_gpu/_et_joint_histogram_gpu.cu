/*
 *  _et_joint_histogram_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_joint_histogram_gpu.h"
#include "_et_joint_histogram_gpu_kernels.cu"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>

void et_joint_histogram_gpu(float **d_array_A, float **d_array_B, int **d_jont_hist, int array_size, int hist_size, float min_A, float max_A, float min_B, float max_B)
{
	//CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_backprojection_size,&backprojection_size,sizeof(int3)));
	
	const unsigned int grid = (unsigned int)ceil(array_size/(float)BLOCK);
	dim3 B(BLOCK,1,1);
	dim3 G(grid,1,1);
	
	et_joint_histogram_gpu_kernel <<<G,B>>> (*d_array_A, *d_array_B, *d_jont_hist, array_size, hist_size, min_A, max_A, min_B, max_B);
	
	CUDA_SAFE_CALL(cudaThreadSynchronize());
}



