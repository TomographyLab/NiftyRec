/*
 *  _tt_line_backproject_ray_gpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _TT_BACKPROJECT_RAY_GPU_H
#define _TT_BACKPROJECT_RAY_GPU_H

#include <_tt_common.h>
// Utilities and System includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <cutil_inline.h>
//#include <cutil_math.h>
#include "_reg_blocksize_gpu.h"
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
//#include <sys/time.h>


extern "C" int set_inViewMatrix(float *invViewMatrix, float_2 detector_scale, float_3 detector_transl, float_3 detector_rotat);
extern "C" void tt_line_backproject_ray_gpu(dim3 gridSize, dim3 blockSize, float *d_projection, float *d_output, uint2 detectorPixels, float3 sourcePosition, uint3 volumeVoxels, float3 volumeSize, float t_step, int interpolation);
extern "C" void copyInvViewMatrix_bk(float *invViewMatrix, size_t sizeofMatrix);

#endif

