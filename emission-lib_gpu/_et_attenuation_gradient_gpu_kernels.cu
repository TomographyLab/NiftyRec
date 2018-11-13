/*
 *  _et_attenuation_gradient_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_attenuation_gradient_gpu.h"

__device__ __constant__ int3 c_backprojection_size;

__global__ void et_attenuation_gradient_gpu_kernel(float *g_activity, float *g_sinogram, float *g_backprojection, float *g_attenuation)
{
        __shared__ float s_sino[BLOCK];

        const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
        const unsigned int pixelNumber = c_backprojection_size.x * c_backprojection_size.y;
        if(tid<pixelNumber){
                //load sinogram to shared mem
                s_sino[threadIdx.x] = g_sinogram[tid];
                // first pass: compute integral attenuation from zero to voxel
                unsigned int index = tid;
                float temp_integral = 0.0f;
                for(unsigned int z=0; z < c_backprojection_size.z; z++){
                        temp_integral += g_attenuation[index];
                        g_backprojection[index] = temp_integral;
                        index += pixelNumber;
                }
                // second pass: compute gradient (integral attenuated activity from horizon to voxel)
                index = tid + pixelNumber*(c_backprojection_size.z-1);
                temp_integral = 0.0f;
                for(unsigned int z=0; z < c_backprojection_size.z; z++){
                        temp_integral += g_activity[index] * exp(-g_backprojection[index]); 
                        g_backprojection[index] = s_sino[threadIdx.x]*temp_integral;
                        index -= pixelNumber;
                }
        }
        return;  
}





