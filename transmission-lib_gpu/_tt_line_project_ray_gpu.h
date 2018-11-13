/*
 *  _tt_line_project_ray_gpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _TT_PROJECT_RAY_GPU_H
#define _TT_PROJECT_RAY_GPU_H

//#include <cutil_inline.h>
#include "_reg_blocksize_gpu.h"
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
//#include <sys/time.h>
#include <_tt_common.h>

extern "C" void setTextureFilterMode(bool bLinearFilter);
extern "C" void initCuda(void *h_volume, cudaExtent volumeSize);
extern "C" void freeCudaBuffers();
extern "C" void tt_line_project_ray_gpu(dim3 gridSize, dim3 blockSize, float *d_output, float3 source_position, float3 volume_size, u_int imageW, u_int imageH, float t_step);
extern "C" void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix);

#endif

