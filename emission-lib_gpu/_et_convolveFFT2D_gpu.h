/*
 *  _et_convolveFFT2D_gpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cufft.h>
//#include <cutil_inline.h>
#include "_reg_blocksize_gpu.h"
#include "nifti1_io.h"


#ifdef __CUDACC__
    typedef float2 fComplex;
#else
    typedef struct{
        float x;
        float y;
    } fComplex;
#endif


int calculateFFTsize(int dataSize);
int et_convolveFFT2D_gpu(float **d_data, int *data_size, float **d_kernel, int *kernel_size, float **d_result);



////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////
//Round a / b to nearest higher integer value
inline int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Align a to nearest higher multiple of b
inline int iAlignUp(int a, int b){
    return (a % b != 0) ?  (a - a % b + b) : a;
}

extern "C" void convolutionClampToBorderCPU(
    float *h_Result,
    float *h_Data,
    float *h_Kernel,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelX,
    int kernelY
);

extern "C" void padKernel(
    float *d_PaddedKernel,
    float *d_Kernel,
    int fftH,
    int fftW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void padDataClampToBorder(
    float *d_PaddedData,
    float *d_Data,
    int fftH,
    int fftW,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void modulateAndNormalize(
    fComplex *d_Dst,
    fComplex *d_Src,
    int fftH,
    int fftW
);

extern "C" void crop_image(
    float *d_Dst,
    float *d_Src,
    int fftH,
    int fftW,
    int imageH,
    int imageW,
    int kernelH,
    int kernelW
);


