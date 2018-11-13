/*
 *  _tt_line_project_ray_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _TTPROJECTRAY_KERNEL_CU_
#define _TTPROJECTRAY_KERNEL_CU_

//#include <cutil_inline.h>
#include "_reg_blocksize_gpu.h"
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <cutil_math.h>
#include <_tt_common.h>

#define MAX_STEPS 1000000000

cudaArray *d_volumeArray = 0;

//texture<VolumeType, 3, cudaReadModeNormalizedFloat> tex;         // 3D texture
texture<VolumeType, 3, cudaReadModeElementType> tex;         // 3D texture

typedef struct {
        float4 m[3];
} float3x4;

__constant__ float3x4 c_invViewMatrix;  // inverse view matrix

struct Ray {
	float3 o;	// origin
	float3 d;	// direction
};

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

__device__
int intersectBox(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invR = make_float3(1.0f) / r.d;
    float3 tbot = invR * (boxmin - r.o);
    float3 ttop = invR * (boxmax - r.o);

    // re-order intersections to find smallest and largest on each axis
    float3 tmin = fminf(ttop, tbot);
    float3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;

    return smallest_tmax > largest_tmin;
}

// transform vector by matrix (no translation)
__device__
float3 mul(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}

// transform vector by matrix with translation
__device__
float4 mul(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}


__global__ void
d_tt_line_project_ray_gpu(float *d_output, float3 sourcePosition, float3 volumeSize, uint imageWidthPixels, uint imageHeightPixels, float tStep, int interpolation)
{
    const uint   image_width_pixels = imageWidthPixels; 
    const uint   image_height_pixels = imageHeightPixels; 
    const float3 volume_size = volumeSize; 
    const float3 source_position = sourcePosition;

    const float  tstep    = tStep;
    const int    maxSteps = MAX_STEPS;          
    //const int    maxSteps = MAX_STEPS;          //(volume_size.x^2+volume_size.y^2+volume_size.z^2)^0.5f/tStep;       //diagonal of the bounding box
    const float3 boxMin = make_float3(0.0f, 0.0f, 0.0f);
    const float3 boxMax = make_float3(volume_size.x, volume_size.y, volume_size.z);

    const float3 rec_volume_size = 1.0f/volume_size;

    //x and y index detector pixel
    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;
    if ((x >= image_width_pixels) || (y >= image_height_pixels)) return;

    //u and v are in normalized detector pixel [0,0]->[1,1]
    float u = (x / (float) image_width_pixels);
    float v = (y / (float) image_height_pixels);

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = source_position;
    //transform and normalize direction vector
    eyeRay.d = normalize(make_float3(mul(c_invViewMatrix, make_float4(u,v,0.0f,1.0f)))-eyeRay.o); 
    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
    if (!hit) return;
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating
    float  sum; 
    float  t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;

    for(int i=0; i<maxSteps; i++) {
        // read from 3D texture
        // remap position to [0, 1] coordinates
        //float sample = tex3D(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
        float sample = tex3D(tex, pos.x*rec_volume_size.x, pos.y*rec_volume_size.y, pos.z*rec_volume_size.z);
        sum = sum + sample;

        t += tstep;
        if (t > tfar) break;

        pos += step;
    }

    d_output[y*image_width_pixels + x] = sum;
}


extern "C" void setTextureFilterMode(bool bLinearFilter)
{
    tex.filterMode = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;
}


extern "C" void initCuda(void *h_volume, cudaExtent volumeSize)
{
    // create 3D array
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    CUDA_SAFE_CALL( cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize) );

    // copy data to 3D array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    CUDA_SAFE_CALL( cudaMemcpy3D(&copyParams) );  

    // set texture parameters
    tex.normalized = true;                      // access with normalized texture coordinates
    tex.filterMode = cudaFilterModeLinear;      // linear interpolation
    tex.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    tex.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    CUDA_SAFE_CALL(cudaBindTextureToArray(tex, d_volumeArray, channelDesc));
}



extern "C" void freeCudaBuffers()
{
    CUDA_SAFE_CALL(cudaFreeArray(d_volumeArray));
}



extern "C" void tt_line_project_ray_gpu(dim3 gridSize, dim3 blockSize, float *d_output, float3 source_position, float3 volume_size, uint imageW, uint imageH, float t_step, int interpolation)
{
	d_tt_line_project_ray_gpu<<<gridSize, blockSize>>>( d_output, source_position, volume_size, imageW, imageH, t_step, interpolation);
}



extern "C" void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
/*    fprintf(stderr,"\nMatrix:");
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[0],invViewMatrix[1],invViewMatrix[3],invViewMatrix[3]);
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[4],invViewMatrix[5],invViewMatrix[6],invViewMatrix[7]);
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[8],invViewMatrix[9],invViewMatrix[10],invViewMatrix[11]);
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[12],invViewMatrix[13],invViewMatrix[14],invViewMatrix[15]);*/
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix) );
}


#endif // #ifndef _TTPROJECTRAY_KERNEL_CU_
