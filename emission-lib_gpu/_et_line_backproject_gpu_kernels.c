/*
 *  _et_line_backproject_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_line_backproject_gpu.h"

__device__ __constant__ int3 c_backprojection_size;


__global__ void et_line_backproject_gpu_kernel(float *g_sinogram, float *g_backprojection)
{
	__shared__ float s_sino[BLOCK];
	
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_backprojection_size.x * c_backprojection_size.y;
	if(tid<pixelNumber){
		//load sinogram to shared mem
		s_sino[threadIdx.x] = g_sinogram[tid];
		
		//backproject
		unsigned int index = tid;
		for(unsigned int z=0; z < c_backprojection_size.z; z++){
			g_backprojection[index] = s_sino[threadIdx.x];
			index += pixelNumber;
		}
	}
	return;  
}


