
__device__ __constant__ int3 c_ImageSize;
__device__ __constant__ float3 c_ImageOrigin;
__device__ __constant__ float3 c_DetectorOrigin;
__device__ __constant__ float2 c_DetectorSize;


__global__ void tt_perspective_positionField_gpu_kernel(float4 *positionFieldArray)
{
	const int tid_x= blockIdx.x*blockDim.x + threadIdx.x;
	const int tid_y= blockIdx.y*blockDim.y + threadIdx.y;

	if(tid_x<c_ImageSize.x)
            {
            if(tid_y<c_ImageSize.z)
                {
		int3 imageSize = c_ImageSize;
		float3 imageOrigin = c_ImageOrigin;
		float3 detectorOrigin = c_DetectorOrigin;
		float2 detectorSize = c_DetectorSize;
		short3 voxelIndex;
                int out_index = imageSize.x * tid_y + tid_x;  
                //int out_index = imageSize.z * tid_x + tid_y; 

		/* The transformation is applied */
		float4 position;
                for(int z_index=0; z_index<imageSize.y; z_index++)
                    {
                    float z_prime = tid_x + imageOrigin.z;
                    float z_ratio = z_prime / detectorOrigin.z;
                    float x_prime = tid_x + detectorOrigin.x * z_ratio;
                    float y_prime = tid_y + detectorOrigin.y * z_ratio;
                    position.x = x_prime; 
                    position.y = y_prime; 
                    position.z = z_prime;
                    position.w = 0.0f;
                    /* the deformation field (real coordinates) is stored */
                    positionFieldArray[out_index] = position;
                    }
                }
	}
}
