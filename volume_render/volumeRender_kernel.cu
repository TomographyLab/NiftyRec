
/*
 * This source file is based on a modification of a file volumeRender_kernel.cu
 * distributed with the NVidia CUDA Toolkit as of 15 Oct. 2013. The original 
 * copyright notice is maintained.  
 *
 * Stefano Pedemonte, Helsinki, Oct. 2013
 */
 
/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */


#ifndef _VOLUMERENDER_KERNEL_CU_
#define _VOLUMERENDER_KERNEL_CU_

#define MAX_STEPS 50000

#include <helper_cuda.h>
#include <helper_math.h>
#include <volumeRender_kernel.h>

static cudaArray *d_volumeArray1 = 0;
static cudaArray *d_volumeArray2 = 0;
static cudaArray *d_transferFuncArray1;
static cudaArray *d_transferFuncArray2;
static cudaArray *d_transferFuncArray1_swap;
static cudaArray *d_transferFuncArray2_swap;

__device__ __constant__ float t_step;
__device__ __constant__ float x_camera;
__device__ __constant__ float y_camera;
__device__ __constant__ float z_camera;

static texture<VolumeType, 3, cudaReadModeNormalizedFloat> tex1;         // 3D texture
static texture<float4, 1, cudaReadModeElementType>         transferTex1; // 1D transfer function texture
unsigned int colormap_inuse1=0;

static texture<VolumeType, 3, cudaReadModeNormalizedFloat> tex2;         // 3D texture
static texture<float4, 1, cudaReadModeElementType>         transferTex2; // 1D transfer function texture
unsigned int colormap_inuse2=0;

//#pragma data_seg (".myseg")
//int temp=0;
//#pragma data_seg()


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

__device__ uint rgbaFloatToInt(float4 rgba)
{
    rgba.x = __saturatef(rgba.x);   // clamp to [0.0, 1.0]
    rgba.y = __saturatef(rgba.y);
    rgba.z = __saturatef(rgba.z);
    rgba.w = __saturatef(rgba.w);
    return (uint(rgba.w*255)<<24) | (uint(rgba.z*255)<<16) | (uint(rgba.y*255)<<8) | uint(rgba.x*255);
}

__global__ void
d_render(uint *d_output, uint imageW, uint imageH,
         float density1, float brightness1, float transferOffset1, float transferScale1, 
         float density2, float brightness2, float transferOffset2, float transferScale2)
{
    const int maxSteps = MAX_STEPS;
    //const float t_step = 0.01f;
    const float opacityThreshold = 0.95f;
    const float3 boxMin = make_float3(-1.0f, -1.0f, -1.0f);
    const float3 boxMax = make_float3(1.0f, 1.0f, 1.0f);

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;
    if ((x >= imageW) || (y >= imageH)) return;

    float u = (x / (float) imageW)*2.0f-1.0f;
    float v = (y / (float) imageH)*2.0f-1.0f;

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(x_camera, y_camera, z_camera, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -(z_camera+2.f)));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);

    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);
    if (!hit) return;
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum1 = make_float4(0.0f);
    float4 sum2 = make_float4(0.0f);
    float4 sum = make_float4(0.0f);
    float  t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*t_step;
    float px, py, pz;
    float sample1, sample2;
    float4 col1, col2;

    float4 col1_prev=make_float4(0.0f);
    float4 col2_prev=make_float4(0.0f);
    float4 col1_prev_prev=make_float4(0.0f);
    float4 col2_prev_prev=make_float4(0.0f);

    for(int i=0; i<maxSteps; i++) {
        // read from 3D texture
        // remap position to [0, 1] coordinates
        px=pos.x*0.5f+0.5f;
        py=pos.y*0.5f+0.5f;
        pz=pos.z*0.5f+0.5f;
        sample1 = tex3D(tex1, px, py, pz);
        sample2 = tex3D(tex2, px, py, pz);

        // lookup in transfer function texture
        col1 = tex1D(transferTex1, (sample1-transferOffset1)*transferScale1);
        col2 = tex1D(transferTex2, (sample2-transferOffset2)*transferScale2);
        col1.w *= density1;
        col2.w *= density2;

        // "under" operator for back-to-front blending
        //sum = lerp(sum, col, col.w);

        // pre-multiply alpha
        col1.x *= col1.w;
        col1.y *= col1.w;
        col1.z *= col1.w;
        col2.x *= col2.w;
        col2.y *= col2.w;
        col2.z *= col2.w;
        // "over" operator for front-to-back blending
        sum1 = sum1 + (0.500*col1+0.333*col1_prev+0.166*col1_prev_prev)*(1.0f - sum1.w);
        sum2 = sum2 + (0.500*col2+0.333*col2_prev+0.166*col2_prev_prev)*(1.0f - sum2.w);

        col1_prev_prev = col1_prev;
        col2_prev_prev = col2_prev;

        col1_prev = col1;
        col2_prev = col2;

        // exit early if opaque
        if (sum1.w > opacityThreshold)
            break;
        if (sum2.w > opacityThreshold)
            break;

        t += t_step;
        if (t > tfar) break;

        pos += step;
    }
    sum1 *= brightness1;
    sum2 *= brightness2;

   sum = sum1+sum2;
    
    // write output color
    d_output[y*imageW + x] = rgbaFloatToInt(sum);
}

extern "C"
void setTextureFilterMode(bool bLinearFilter)
{
    tex1.filterMode = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;
    tex2.filterMode = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;
}

extern "C"
void initCuda(uint w, uint h, uint d)
{
    cudaExtent volumeSize = make_cudaExtent(w,h,d);
    // create 3D array
    cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<VolumeType>();
    cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<VolumeType>();
    cudaMalloc3DArray(&d_volumeArray1, &channelDesc1, volumeSize);
    cudaMalloc3DArray(&d_volumeArray2, &channelDesc2, volumeSize);

    // set texture parameters
    tex1.normalized = true;                      // access with normalized texture coordinates
    tex1.filterMode = cudaFilterModeLinear;      // linear interpolation
    tex1.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    tex1.addressMode[1] = cudaAddressModeClamp;

    tex2.normalized = true;                      // access with normalized texture coordinates
    tex2.filterMode = cudaFilterModeLinear;      // linear interpolation
    tex2.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    tex2.addressMode[1] = cudaAddressModeClamp;


    // bind array to 3D texture
    cudaBindTextureToArray(tex1, d_volumeArray1, channelDesc1);
    cudaBindTextureToArray(tex2, d_volumeArray2, channelDesc2);

    // create transfer function texture
    float4 transferFunc1[] = {
        {  0.0, 0.0, 0.0, 0.0, },
        {  0.6, 0.6, 0.6, 0.2, },
        {  0.0, 1.0, 0.0, 0.5, },
        {  1.0, 0.0, 0.0, 0.5, }, };

    float4 transferFunc2[] = {
        {  0.0, 0.0, 0.0, 0.0, },
        {  0.2, 0.2, 0.2, 0.6, },
        {  0.4, 0.4, 0.4, 0.6, },
        {  0.6, 0.6, 0.6, 0.6, },
        {  0.8, 0.8, 0.8, 0.6, },
        {  1.0, 1.0, 1.0, 0.6, },
    };

    cudaChannelFormatDesc channelDesc1a = cudaCreateChannelDesc<float4>();
    cudaArray* d_transferFuncArray1;
    cudaMallocArray( &d_transferFuncArray1, &channelDesc1a, sizeof(transferFunc1)/sizeof(float4), 1); 
    cudaMemcpyToArray( d_transferFuncArray1, 0, 0, transferFunc1, sizeof(transferFunc1), cudaMemcpyHostToDevice);
    transferTex1.filterMode = cudaFilterModeLinear;
    transferTex1.normalized = true;    // access with normalized texture coordinates
    transferTex1.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    cudaChannelFormatDesc channelDesc2a = cudaCreateChannelDesc<float4>();
    cudaArray* d_transferFuncArray2;
    cudaMallocArray( &d_transferFuncArray2, &channelDesc2a, sizeof(transferFunc2)/sizeof(float4), 1); 
    cudaMemcpyToArray( d_transferFuncArray2, 0, 0, transferFunc2, sizeof(transferFunc2), cudaMemcpyHostToDevice);
    transferTex2.filterMode = cudaFilterModeLinear;
    transferTex2.normalized = true;    // access with normalized texture coordinates
    transferTex2.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    cudaBindTextureToArray( transferTex1, d_transferFuncArray1, channelDesc1a);
    cudaBindTextureToArray( transferTex2, d_transferFuncArray2, channelDesc2a);
}

extern "C" 
void freeCudaBuffers()
{
    cudaFreeArray(d_volumeArray1);
    if(colormap_inuse1==0) {
        cudaFreeArray(d_transferFuncArray1);}
    else {
        cudaFreeArray(d_transferFuncArray1_swap);
        cudaFreeArray(d_volumeArray2); }
    if(colormap_inuse2==0) {
        cudaFreeArray(d_transferFuncArray2); }
    else {
        cudaFreeArray(d_transferFuncArray2_swap); }
}


extern "C"
void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint imageW, uint imageH, 
				   float density1, float brightness1, float transferOffset1, float transferScale1, 
                                   float density2, float brightness2, float transferOffset2, float transferScale2)
{
	d_render<<<gridSize, blockSize>>>( d_output, imageW, imageH, 
                                   density1, brightness1, transferOffset1, transferScale1, 
                                   density2, brightness2, transferOffset2, transferScale2);
}

extern "C"
void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
    cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix);
}

extern "C" void setTexture1(VolumeType* h_volume,uint w,uint h,uint d)
{
    cudaExtent volumeSize = make_cudaExtent(w,h,d);
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray1;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);  
}

extern "C"
void setTexture2(VolumeType* h_volume,uint w,uint h,uint d)
{
    cudaExtent volumeSize = make_cudaExtent(w,h,d);
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray2;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);  
}


extern "C"
void setTransferFunction1(float *colormap, unsigned int colormap_size)
{
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
    if (colormap_inuse1==0) {
        colormap_inuse1=1;
        //allocate new colormap device array and copy colormap to device
        cudaMallocArray( &d_transferFuncArray1_swap, &channelDesc, colormap_size, 1); 
        cudaMemcpyToArray( d_transferFuncArray1_swap, 0, 0, colormap, colormap_size*4*sizeof(float), cudaMemcpyHostToDevice);    
        //unbind previous colormap and bind texture to the new one
        cudaUnbindTexture( transferTex1 );
        transferTex1.filterMode = cudaFilterModeLinear;
        transferTex1.normalized = true;    // access with normalized texture coordinates
        transferTex1.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
        cudaBindTextureToArray( transferTex1, d_transferFuncArray1_swap, channelDesc);
        cudaFreeArray(d_transferFuncArray1);
        }
    else {
        colormap_inuse1=0;
        //allocate new colormap device array and copy colormap to device
        cudaMallocArray( &d_transferFuncArray1, &channelDesc, colormap_size, 1); 
        cudaMemcpyToArray( d_transferFuncArray1, 0, 0, colormap, colormap_size*4*sizeof(float), cudaMemcpyHostToDevice);   
        //unbind previous colormap and bind texture to the new one
        cudaUnbindTexture( transferTex1 );
        transferTex1.filterMode = cudaFilterModeLinear;
        transferTex1.normalized = true;    // access with normalized texture coordinates
        transferTex1.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
        cudaBindTextureToArray( transferTex1, d_transferFuncArray1, channelDesc);
        //free previous colormap
        cudaFreeArray(d_transferFuncArray1_swap);
        } 
}

extern "C"
void setTransferFunction2(float *colormap, unsigned int colormap_size)
{
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
    if (colormap_inuse2==0) {
        colormap_inuse2=1;
        //allocate new colormap device array and copy colormap to device
        cudaMallocArray( &d_transferFuncArray2_swap, &channelDesc, colormap_size, 1); 
        cudaMemcpyToArray( d_transferFuncArray2_swap, 0, 0, colormap, colormap_size*4*sizeof(float), cudaMemcpyHostToDevice);    
        //unbind previous colormap and bind texture to the new one
        cudaUnbindTexture( transferTex2 );
        transferTex2.filterMode = cudaFilterModeLinear;
        transferTex2.normalized = true;    // access with normalized texture coordinates
        transferTex2.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
        cudaBindTextureToArray( transferTex2, d_transferFuncArray2_swap, channelDesc);
        cudaFreeArray(d_transferFuncArray2);
        }
    else {
        colormap_inuse2=0;
        //allocate new colormap device array and copy colormap to device
        cudaMallocArray( &d_transferFuncArray2, &channelDesc, colormap_size, 1); 
        cudaMemcpyToArray( d_transferFuncArray2, 0, 0, colormap, colormap_size*4*sizeof(float), cudaMemcpyHostToDevice);   
        //unbind previous colormap and bind texture to the new one
        cudaUnbindTexture( transferTex2 );
        transferTex1.filterMode = cudaFilterModeLinear;
        transferTex1.normalized = true;    // access with normalized texture coordinates
        transferTex1.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
        cudaBindTextureToArray( transferTex2, d_transferFuncArray2, channelDesc);
        //free previous colormap
        cudaFreeArray(d_transferFuncArray2_swap);
        } 

}

extern "C"
void setRayStep(float step)
{
    cudaMemcpyToSymbol(t_step, &step, sizeof(float));
}

extern "C"
void setCameraPosition(float x, float y, float z)
{
    cudaMemcpyToSymbol(x_camera, &x, sizeof(float));
    cudaMemcpyToSymbol(y_camera, &y, sizeof(float));
    cudaMemcpyToSymbol(z_camera, &z, sizeof(float));
}

#endif // #ifndef _VOLUMERENDER_KERNEL_CU_
