/*
 *  _et_line_integrals_attenuated_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_line_integral_attenuated_gpu.h"

__device__ __constant__ int3 c_ImageSize;


__global__ void et_line_integral_attenuated_gpu_kernel(float *g_activity, float *g_attenuation, float *g_background_image, float *g_sinogram, float *g_partialsum, float background_activity) 
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.y;
	if(tid<pixelNumber){
		unsigned int index=tid;
                unsigned int index_inv=tid;
                float sum_attenuation=0.0f;
		float sum_activity=0.0f;
                if (g_partialsum==NULL)
                    {
                    if (g_activity==NULL && g_attenuation!=NULL)
                        {
                        for(unsigned int z=0; z<c_ImageSize.z; z++)
                            {
                            sum_attenuation += g_attenuation[index];
                            index += pixelNumber;
                            }
                        }
                    else if (g_activity!=NULL && g_attenuation!=NULL)
                        {
		        for(unsigned int z=0; z<c_ImageSize.z; z++)
                            {
                            sum_attenuation += g_attenuation[index];
                            sum_activity    += g_activity[index]*exp(-sum_attenuation);
                            index += pixelNumber;
                            }
                        }
                    else if (g_activity!=NULL && g_attenuation==NULL)
                        {
                        for(unsigned int z=0; z<c_ImageSize.z; z++)
                            {
                            sum_activity += g_activity[index];
                            index += pixelNumber;
                            }
                        }
                    else
                        return;
                    sum_activity += background_activity*exp(-sum_attenuation);
                    g_sinogram[tid]=sum_activity;
                    }
                else
                    {
                    if (g_activity==NULL && g_attenuation!=NULL) //OK
                        {
                        index_inv = tid+pixelNumber*(c_ImageSize.z-1);
                        sum_activity = background_activity; 
                        if (g_background_image!=NULL)
                            sum_activity = sum_activity*g_background_image[tid];
		        for(unsigned int z=0; z<c_ImageSize.z; z++)
                            {
                            g_partialsum[index_inv] = sum_activity; 
                            sum_activity = sum_activity * exp(-g_attenuation[index_inv]); 
                            index_inv -= pixelNumber;
                            }
                        g_sinogram[tid]=sum_activity;
                        }
                    else if (g_activity!=NULL && g_attenuation!=NULL) // OK
                        {
                        index_inv = tid+pixelNumber*(c_ImageSize.z-1);
                        sum_activity = background_activity; 
                        if (g_background_image!=NULL)
                            sum_activity = sum_activity*g_background_image[tid];
		        for(unsigned int z=0; z<c_ImageSize.z; z++)
                            {
                            g_partialsum[index_inv] = sum_activity; 
                            sum_activity = (sum_activity+g_activity[index_inv]) * exp(-g_attenuation[index_inv]); 
                            index_inv -= pixelNumber;
                            }
                        g_sinogram[tid]=sum_activity;
                        }
                    else if (g_activity!=NULL && g_attenuation==NULL) // OK
                        {
                        index_inv = tid+pixelNumber*(c_ImageSize.z-1);
                        sum_activity = background_activity; 
                        if (g_background_image!=NULL)
                            sum_activity = sum_activity*g_background_image[tid];
		        for(unsigned int z=0; z<c_ImageSize.z; z++)
                            {
                            g_partialsum[index_inv] = sum_activity; 
                            sum_activity = (sum_activity+g_activity[index_inv]); 
                            index_inv -= pixelNumber;
                            }
                        g_sinogram[tid]=sum_activity;
                        }
                    else
                        return;
                    }
	}
	return; 	
}




