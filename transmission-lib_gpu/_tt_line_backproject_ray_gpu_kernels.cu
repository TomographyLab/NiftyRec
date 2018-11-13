/*
 *  _tt_line_backproject_ray_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _TTBACKPROJECTRAY_KERNEL_CU_
#define _TTBACKPROJECTRAY_KERNEL_CU_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <cutil_inline.h>
#include <cutil_math.h>
#include "device_functions.h"
#include <_tt_line_backproject_ray_gpu.h>

#define MAX_STEPS 1000000000

typedef struct {
        float4 m[3];
} float3x4;

__constant__ float3x4 c_invViewMatrix_bk;  // inverse view matrix

struct Ray {
	float3 o;	// origin
	float3 d;	// direction
};

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

__device__
int intersectBox_(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
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
float3 mul_(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}

// transform vector by matrix with translation
__device__
float4 mul_(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}


__global__ void
d_tt_line_backproject_ray_gpu(float *d_projection, float *d_output, uint3 volumeVoxels, float3 sourcePosition, float3 volumeSize, uint2 detectorPixels, float tStep, int interpolation)
{
    const uint   image_width_pixels  = detectorPixels.x; 
    const uint   image_height_pixels = detectorPixels.y; 
    const float3 volume_size = volumeSize; 
    const uint3  volume_voxels = volumeVoxels;
    const float3 source_position = sourcePosition;

    const float  tstep    = tStep;
    const int    maxSteps = MAX_STEPS; //(volume_size.x^2+volume_size.y^2+volume_size.z^2)^0.5f/tStep;       //diagonal of the bounding box
    const float3 boxMin = make_float3(0.0f, 0.0f, 0.0f);
    const float3 boxMax = make_float3(volume_size.x, volume_size.y, volume_size.z);

    //x and y index detector pixel
    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;
    if ((x >= image_width_pixels) || (y >= image_height_pixels)) return;
//    if(!(x%2==0 && y%2==0)) return;
    //u and v are in normalized detector pixel [0,0]->[1,1]
    float u = (x / (float) image_width_pixels);
    float v = (y / (float) image_height_pixels);

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = source_position;
    //transform and normalize direction vector
    eyeRay.d = normalize(make_float3(mul_(c_invViewMatrix_bk, make_float4(u,v,0.0f,1.0f)))-eyeRay.o); 
    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox_(eyeRay, boxMin, boxMax, &tnear, &tfar);
    if (!hit) return;
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane
    int index; 
    float proj = d_projection[y*image_width_pixels + x];

    // march along ray from front to back
    float  t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;


    for(int ii=0; ii<maxSteps; ii++) 
    {
        if (t > tfar-(volume_size.x/100)) break;
        if (interpolation==INTERP_NEAREST)
        {
            index = (int) floor(pos.z)*volume_voxels.x*volume_voxels.y + floor(pos.y)*volume_voxels.x + floor(pos.x);
    //        atomicAdd(d_output+index, proj);
            if (pos.x<volume_voxels.x && pos.x>=0)
                if (pos.y<volume_voxels.y && pos.y>=0)
                    if (pos.z<volume_voxels.z && pos.z>=0)
            d_output[index] += proj;        
        }
        else if (interpolation==INTERP_TRILINEAR)
        {
            if (pos.x<(volume_voxels.x-1) && pos.x>=1) 
            {
                if (pos.y<(volume_voxels.y-1) && pos.y>=1) 
                {
                    if (pos.z<(volume_voxels.z-1) && pos.z>=1) 
                    {
                        int3 p000, p001, p010, p011, p100, p101, p110, p111;
                        p000.x = (int)floor(pos.x); p000.y = (int)floor(pos.y); p000.z = (int)floor(pos.z); 
                        p001 = p000; p001.x+=1;
                        p010 = p000; p010.y+=1;
                        p011 = p010; p011.x+=1;
                        p100 = p000; p100.z+=1;
                        p101 = p001; p101.z+=1;
                        p110 = p100; p110.y+=1;
                        p111 = p110; p111.x+=1; 
                        float3 d;
                        d.x = pos.x-p000.x;
                        d.y = pos.y-p000.y;
                        d.z = pos.z-p000.z; 
                        float w000, w001, w010, w011, w100, w101, w110, w111;                      
                        w000 = (1-d.z)*(1-d.y)*(1-d.x)*proj;
                        w001 = (1-d.z)*(1-d.y)*( d.x )*proj;
                        w010 = (1-d.z)*( d.y )*(1-d.x)*proj;
                        w011 = (1-d.z)*( d.y )*( d.x )*proj;
                        w100 = ( d.z )*(1-d.y)*(1-d.x)*proj;
                        w101 = ( d.z )*(1-d.y)*( d.x )*proj;
                        w110 = ( d.z )*( d.y )*(1-d.x)*proj;
                        w111 = ( d.z )*( d.y )*( d.x )*proj;
                        d_output[p000.z*volumeVoxels.y*volumeVoxels.x+p000.y*volumeVoxels.x+p000.x] += w000;
                        d_output[p001.z*volumeVoxels.y*volumeVoxels.x+p001.y*volumeVoxels.x+p001.x] += w001;
                        d_output[p010.z*volumeVoxels.y*volumeVoxels.x+p010.y*volumeVoxels.x+p010.x] += w010;
                        d_output[p011.z*volumeVoxels.y*volumeVoxels.x+p011.y*volumeVoxels.x+p011.x] += w011;
                        d_output[p100.z*volumeVoxels.y*volumeVoxels.x+p100.y*volumeVoxels.x+p100.x] += w100;
                        d_output[p101.z*volumeVoxels.y*volumeVoxels.x+p101.y*volumeVoxels.x+p101.x] += w101;
                        d_output[p110.z*volumeVoxels.y*volumeVoxels.x+p110.y*volumeVoxels.x+p110.x] += w110;
                        d_output[p111.z*volumeVoxels.y*volumeVoxels.x+p111.y*volumeVoxels.x+p111.x] += w111;
                    }
                }
            }
        }

        __syncthreads();

        t += tstep;
        pos += step;

    }
}




extern "C" void tt_line_backproject_ray_gpu(dim3 gridSize, dim3 blockSize, float *d_projection, float *d_output, uint2 detectorPixels, float3 sourcePosition, uint3 volumeVoxels, float3 volumeSize, float t_step, int interpolation)
{
//        fprintf(stderr,"\nGrid:  %d %d %d",gridSize.x,gridSize.y,gridSize.z);
//        fprintf(stderr,"\nBlock: %d %d %d",blockSize.x,blockSize.y,blockSize.z);
//        fprintf(stderr,"\nVolume size:    %f %f %f",volumeSize.x,volumeSize.y,volumeSize.z);
//        fprintf(stderr,"\nVolume voxels:  %d %d %d",volumeVoxels.x,volumeVoxels.y,volumeVoxels.z);
//        fprintf(stderr,"\nSource:         %f %f %f",sourcePosition.x,sourcePosition.y,sourcePosition.z);
//        fprintf(stderr,"\ndetectorPixels: %d %d   ",detectorPixels.x, detectorPixels.y);
	d_tt_line_backproject_ray_gpu<<<gridSize, blockSize>>>( d_projection, d_output, volumeVoxels, sourcePosition, volumeSize, detectorPixels, t_step, interpolation);
        CUDA_SAFE_CALL(cudaThreadSynchronize());
}



extern "C" void copyInvViewMatrix_bk(float *invViewMatrix, size_t sizeofMatrix)
{
/*    fprintf(stderr,"\nMatrix:");
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[0],invViewMatrix[1],invViewMatrix[3],invViewMatrix[3]);
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[4],invViewMatrix[5],invViewMatrix[6],invViewMatrix[7]);
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[8],invViewMatrix[9],invViewMatrix[10],invViewMatrix[11]);
    fprintf(stderr,"\n %4.2f %4.2f %4.2f %4.2f ",invViewMatrix[12],invViewMatrix[13],invViewMatrix[14],invViewMatrix[15]); */
    CUDA_SAFE_CALL( cudaMemcpyToSymbol("c_invViewMatrix_bk", invViewMatrix, sizeofMatrix, 0, cudaMemcpyHostToDevice) );
//    CUDA_SAFE_CALL(cudaMemcpy(c_invViewMatrix_bk, invViewMatrix, sizeofMatrix, cudaMemcpyHostToDevice));
}




#endif // #ifndef _TTBACKPROJECTRAY_KERNEL_CU_

