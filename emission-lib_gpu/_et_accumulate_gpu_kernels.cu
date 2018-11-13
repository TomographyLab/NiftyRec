/*
 *  _et_attenuation_gradient_gpu_kernels.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et_accumulate_gpu.h"

__device__ __constant__ int3 c_image_size;

__global__ void et_accumulate_gpu_kernel(float *g_A, float *g_B)
{
	__shared__ float s_A[BLOCK];
	__shared__ float s_B[BLOCK];
	
	const unsigned int index  = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int tid    = threadIdx.x;
	const unsigned int n_pixels = c_image_size.x * c_image_size.y * c_image_size.z;
	if(index < n_pixels){
		s_A[tid] = g_A[index];
		s_B[tid] = g_B[index];
		__syncthreads();
		s_B[tid] = s_A[tid] + s_B[tid];
		__syncthreads();
		g_B[index] = s_B[tid];
	}
	return;  
}


