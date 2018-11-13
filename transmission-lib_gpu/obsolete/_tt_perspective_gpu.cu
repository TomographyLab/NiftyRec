#include "_tt_perspective_gpu.h"
#include "_tt_perspective_gpu_kernels.cu"

#define BLOCK 16

void tt_perspective_positionField_gpu(nifti_image *attenuation, float *image_origin, float *detector_origin, float *detector_size, float4 **d_positionField)
{
	int3 imageSize = make_int3(attenuation->dim[1],attenuation->dim[2],attenuation->dim[3]);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ImageSize,&imageSize,sizeof(int3)));

        float3 imageOrigin = make_float3(image_origin[0],image_origin[1],image_origin[2]);
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ImageOrigin,&imageOrigin,sizeof(float3)));

        float3 detectorOrigin = make_float3(detector_origin[0],detector_origin[1],detector_origin[2]);
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_DetectorOrigin,&detectorOrigin,sizeof(float3)));

        float2 detectorSize = make_float2(detector_size[0],detector_size[1]);
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_DetectorSize,&detectorSize,sizeof(float2)));

	
	//const unsigned int Grid = (unsigned int)ceil(img->nx*img->ny/(float)BLOCK);
	const unsigned int grid_x = (unsigned int)ceil(attenuation->dim[1]/(float)BLOCK);
        const unsigned int grid_y = (unsigned int)ceil(attenuation->dim[3]/(float)BLOCK);
	dim3 B1(BLOCK,BLOCK,1);
	dim3 G1(grid_x,grid_y,1);
	
//	float *currentCamPointer = (*d_sinogram) + cam * img->dim[1] * img->dim[3] ;
	
	tt_perspective_positionField_gpu_kernel <<<G1,B1>>> (*d_positionField);

	CUDA_SAFE_CALL(cudaThreadSynchronize());
}


