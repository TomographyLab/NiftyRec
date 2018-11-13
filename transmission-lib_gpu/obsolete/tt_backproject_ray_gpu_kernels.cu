
#ifndef _TTBACKPROJECTRAY_KERNEL_CU_
#define _TTBACKPROJECTRAY_KERNEL_CU_

#include <cutil_inline.h>
#include <cutil_math.h>
#include <tt_backproject_ray_gpu.h>
#include "device_functions.h"

#define MAX_STEPS 10000

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
d_tt_backproject_ray(float *d_projection, float *d_output, uint3 volumeVoxels, float3 sourcePosition, float3 volumeSize, uint imageWidthPixels, uint imageHeightPixels, float tStep)
{
    const uint   image_width_pixels  = imageWidthPixels; 
    const uint   image_height_pixels = imageHeightPixels; 
    const float3 volume_size = volumeSize; 
    const uint3  volume_voxels = volumeVoxels;
    const float3 source_position = sourcePosition;

    const float  tstep    = tStep;
    const int    maxSteps = MAX_STEPS; //(volume_size.x^2+volume_size.y^2+volume_size.z^2)^0.5f/tStep;       //diagonal of the bounding box
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

    //load projection from global memory
    float proj = d_projection[y*image_width_pixels + x];

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = source_position;
    //transform and normalize direction vector
    eyeRay.d = normalize(make_float3(mul(c_invViewMatrix, make_float4(u,v,0.0f,1.0f)))-eyeRay.o); 
    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
//    if (!hit) return;
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back
    float  t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;

    for(int i=0; i<maxSteps; i++) 
    {
        uint X = (pos.x*rec_volume_size.x)*volume_voxels.x; 
        uint Y = (pos.y*rec_volume_size.y)*volume_voxels.y;
        uint Z = (pos.z*rec_volume_size.z)*volume_voxels.z;
        if (X>=volume_voxels.x-1)
            X = volume_voxels.x-1;
        if (Y>=volume_voxels.y-1)
            Y = volume_voxels.y-1;
        if (Z>=volume_voxels.z-1)
            Z = volume_voxels.z-1;

        
        //atomicAdd(d_output, proj);
        d_output[Z*volume_voxels.x*volume_voxels.y + Y*volume_voxels.x + X] = Z*volume_voxels.x*volume_voxels.y + Y*volume_voxels.x + X;
//d_output[Z*volume_voxels.x*volume_voxels.y + Y*volume_voxels.x + X]+proj;        
        __syncthreads();


        t += tstep;
        if (t > tfar) break;

        pos += step;
    }
}



extern "C"
void tt_backproject_ray_kernel(dim3 gridSize, dim3 blockSize, float *d_projection, float *d_output, uint3 volumeVoxels, float3 source_position, float3 volume_size, uint imageW, uint imageH, float t_step)
{
//        fprintf(stderr,"\nGrid:  %d %d %d",gridSize.x,gridSize.y,gridSize.z);
//        fprintf(stderr,"\nBlock: %d %d %d",blockSize.x,blockSize.y,blockSize.z);
//        fprintf(stderr,"\nVolume size:   %f %f %f",volume_size.x,volume_size.y,volume_size.z);
//        fprintf(stderr,"\nVolume voxels: %d %d %d",volumeVoxels.x,volumeVoxels.y,volumeVoxels.z);
//        fprintf(stderr,"\nSource:        %f %f %f",source_position.x,source_position.y,source_position.z);
	d_tt_backproject_ray<<<gridSize, blockSize>>>( d_projection, d_output, volumeVoxels, source_position, volume_size, imageW, imageH, t_step);
        CUDA_SAFE_CALL(cudaThreadSynchronize());
}


extern "C"
void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
    cutilSafeCall( cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix) );
}




#endif // #ifndef _TTBACKPROJECTRAY_KERNEL_CU_
