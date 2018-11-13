/*
 *  _pet_line_backproject_compressed_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Harvard University, Martinos Center for Biomedical Imaging
 *  Jan. 2014.
 */

#include "_pet_line_backproject_compressed_gpu.h"
#include "_pet_line_backproject_compressed_gpu_kernels.cu"

void pet_line_backproject_compressed_gpu(float *d_backprojection, float *d_attenuation, float *d_projection, unsigned short *d_locations, unsigned int N_locations, unsigned int N_u, unsigned int N_v, unsigned int N_samples, unsigned int direction, unsigned int block_size)
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_locations, &N_locations,sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_u,         &N_u,sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_v,         &N_v,sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_N_samples,   &N_samples,sizeof(unsigned int)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_block_size,  &block_size,sizeof(unsigned int)));

//fprintf(stderr,"### pet_line_backproject_compressed_gpu -  Nu: %d   Nv: %d   N_samples: %d   N_locations: %d   \n",N_u,N_v,N_samples,N_locations);

    const unsigned int Grid = (unsigned int)ceil((float)N_locations/(float)block_size); 
    
    // note: block_size affects the performance, a good value is 512
    
	dim3 B1(block_size,1,1);
	dim3 G1(Grid,1,1);

	pet_line_backproject_compressed_gpu_kernel <<<G1,B1>>> (d_backprojection, d_attenuation, d_projection, d_locations, direction);

	CUDA_SAFE_CALL(cudaThreadSynchronize()); 
}



