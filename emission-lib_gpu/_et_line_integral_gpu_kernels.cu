/*
 *  _et_line_integral_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_line_integral_gpu.h"

__device__ __constant__ int3 c_ImageSize;


__global__ void et_line_integral_gpu_kernel(float *g_activity, float *g_sinogram, float background_activity)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.y;
	if(tid<pixelNumber){
		unsigned int index=tid;
		float sum=0.0f;
		for(unsigned int z=0; z<c_ImageSize.z; z++)
                        {
			sum += g_activity[index];
			index += pixelNumber;
			}
                sum += background_activity;
		g_sinogram[tid]=sum;
	}
	return; 	
}



