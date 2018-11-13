/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "_et_convolveFFT2D_gpu.h"

#if(0)
    texture<float, 1, cudaReadModeElementType> texFloat;
    size_t offset = 0;
    #define   LOAD_FLOAT(i) tex1Dfetch(texFloat, i)
    #define  SET_FLOAT_BASE CUDA_SAFE_CALL( cudaBindTexture(&offset, texFloat, d_Src, fftH*fftW*sizeof(float)) )

    #define  UNSET_FLOAT_BASE cudaUnbindTexture( texFloat );
#else
    #define  LOAD_FLOAT(i) d_Src[i]
    #define SET_FLOAT_BASE

    #define UNSET_FLOAT_BASE
#endif


////////////////////////////////////////////////////////////////////////////////
/// Position convolution kernel center at (0, 0) in the image
////////////////////////////////////////////////////////////////////////////////
__global__ void padKernel_kernel(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
){
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int x = blockDim.x * blockIdx.x + threadIdx.x;

    if(y < kernelH && x < kernelW){
        int ky = y - kernelY; if(ky < 0) ky += fftH;
        int kx = x - kernelX; if(kx < 0) kx += fftW;
        d_Dst[ky * fftW + kx] = LOAD_FLOAT(y * kernelW + x);
    }
}

extern "C" void padKernel(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
){
    assert(d_Src != d_Dst);
    dim3 threads(32, 8);
    dim3 grid(iDivUp(kernelW, threads.x), iDivUp(kernelH, threads.y));

    SET_FLOAT_BASE;
    padKernel_kernel<<<grid, threads>>>(
        d_Dst,
        d_Src,
        fftH,
        fftW,
        kernelH,
        kernelW,
        kernelY,
        kernelX
    );
    CUDA_CHECK_MSG("padKernel_kernel<<<>>> execution failed\n");
UNSET_FLOAT_BASE;
}


////////////////////////////////////////////////////////////////////////////////
// Prepare data for "pad to border" addressing mode
////////////////////////////////////////////////////////////////////////////////
__global__ void padDataClampToBorder_kernel(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
){
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int x = blockDim.x * blockIdx.x + threadIdx.x;
    const int borderH = dataH + kernelY;
    const int borderW = dataW + kernelX;

    if(y < fftH && x < fftW){
        int dy, dx;

        if(y < dataH) dy = y;
        if(x < dataW) dx = x;
        if(y >= dataH && y < borderH) dy = dataH - 1;
        if(x >= dataW && x < borderW) dx = dataW - 1;
        if(y >= borderH) dy = 0;
        if(x >= borderW) dx = 0;

        d_Dst[y * fftW + x] = LOAD_FLOAT(dy * dataW + dx);
    }
}

extern "C" void padDataClampToBorder(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int dataH,
    int dataW,
    int kernelW,
    int kernelH,
    int kernelY,
    int kernelX
){
    assert(d_Src != d_Dst);
    dim3 threads(32, 8);
    dim3 grid(iDivUp(fftW, threads.x), iDivUp(fftH, threads.y));

    SET_FLOAT_BASE;
    padDataClampToBorder_kernel<<<grid, threads>>>(
        d_Dst,
        d_Src,
        fftH,
        fftW,
        dataH,
        dataW,
        kernelH,
        kernelW,
        kernelY,
        kernelX
    );
    CUDA_CHECK_MSG("padDataClampToBorder_kernel<<<>>> execution failed\n");
UNSET_FLOAT_BASE;
}


////////////////////////////////////////////////////////////////////////////////
// Modulate Fourier image of padded data by Fourier image of padded kernel
// and normalize by FFT size
////////////////////////////////////////////////////////////////////////////////
__global__ void modulateAndNormalizeKernel(
    fComplex *d_Dst,
    fComplex *d_Src,
    int dataSize,
    float c
){
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i >= dataSize)
        return;

    fComplex a = d_Src[i];
    fComplex b = d_Dst[i];
    fComplex R = {c * (a.x * b.x - a.y * b.y), c * (a.y * b.x + a.x * b.y)};

    d_Dst[i] = R;
}

extern "C" void modulateAndNormalize(
    fComplex *d_Dst,
    fComplex *d_Src,
    int fftH,
    int fftW
){
    int dataSize = fftH * (fftW / 2 + 1);
    modulateAndNormalizeKernel<<<iDivUp(dataSize, 256), 256>>>(
        d_Dst,
        d_Src,
        dataSize,
        1.0f / (float)(fftW * fftH)
    );
    CUDA_CHECK_MSG("modulateAndNormalize() execution failed\n");
}



////////////////////////////////////////////////////////////////////////////////
// Crop image 
////////////////////////////////////////////////////////////////////////////////
__global__ void crop_image_kernel(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int imageH,
    int imageW,
    int kernelH,
    int kernelW
){
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int x = blockDim.x * blockIdx.x + threadIdx.x;

    if(y < imageH && x < imageW){
        d_Dst[y*imageW+x] = LOAD_FLOAT(y*(fftW)+x); //FIXME
        }
}

extern "C" void crop_image(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int imageH,
    int imageW,
    int kernelH,
    int kernelW
){
    assert(d_Src != d_Dst);
    dim3 threads(32, 8);
    dim3 grid(iDivUp(imageW, threads.x), iDivUp(imageH, threads.y));

    SET_FLOAT_BASE;
    crop_image_kernel<<<grid, threads>>>(
        d_Dst,
        d_Src,
        fftH,
        fftW,
        imageH,
        imageW,
        kernelH,
        kernelW
    );
    CUDA_CHECK_MSG("padKernel_kernel<<<>>> execution failed\n");
UNSET_FLOAT_BASE;
}


