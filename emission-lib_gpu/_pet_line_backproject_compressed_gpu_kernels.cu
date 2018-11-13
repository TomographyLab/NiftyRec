/*
 *  _pet_line_backproject_compressed_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Harvard University, Martinos Center for Biomedical Imaging
 *  Jan. 2014.
 */

#include "_pet_line_backproject_compressed_gpu.h"

__device__ __constant__ int  c_N_u; 
__device__ __constant__ int  c_N_v;
__device__ __constant__ int  c_N_locations; 
__device__ __constant__ int  c_N_samples;
__device__ __constant__ int  c_block_size;

__global__ void pet_line_backproject_compressed_gpu_kernel(float *g_backprojection, float *g_attenuation, float *g_projection, unsigned short *g_locations, unsigned int direction)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int N_locations = c_N_locations; 
	const unsigned int N_samples   = c_N_samples; 
	const unsigned int N_u         = c_N_u; 
	const unsigned int N_v         = c_N_v; 
	const unsigned int pixelNumber = N_u * N_v;

//    __shared__ float s_projection[512];	

    if(tid<N_locations){ 
        unsigned short u   = g_locations[tid*3];
        unsigned short v   = g_locations[tid*3+1]; 
        unsigned int index; 
		//load projection to shared mem
//		s_projection[threadIdx.x] = g_projection[tid];

        if (direction==1)
            index = u  +  v*N_u; 
        else if (direction==2)
            index = u*N_v  +  v; 
        else if (direction==3)
            index = u*N_samples  +  v*N_u*N_samples;
        else if (direction==4)
            index = u*N_samples  +  v*N_samples*N_u;
        else if (direction==5)
            index = u*N_v*N_samples  +  v;
        else if (direction==6)
            index = u  +  v*N_u*N_samples;
        else 
            index = u*N_v  +  v; 
 
        for(unsigned int z=0; z<N_samples; z++) {     
//              g_backprojection[index] = s_projection[threadIdx.x];
              g_backprojection[index] = g_projection[tid];
       
            if (direction==1)
                index += pixelNumber; 
            else if (direction==2)
                index += pixelNumber;  
            else if (direction==3)
                index += 1;
            else if (direction==4)
                index += 1;
            else if (direction==5)
                index += N_v;              
            else if (direction==6)
                index += N_u;   
            else 
                index += pixelNumber;  
        }
	}
	return; 
}	



