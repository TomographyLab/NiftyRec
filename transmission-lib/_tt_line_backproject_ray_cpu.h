/*
 *  _tt_line_backproject_ray_cpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, 2019-2013.
 *  CMIC - Centre for Medical Image Computing, UCL, London. (2009-2013)
 *  Aalto University School of Science, Helsinki. (2013) 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <_tt_common.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//FIXME complete the CPU version (no CUDA)
#include <_tt_common.h>

extern "C" int tt_line_backproject_ray_cpu(float *out_backprojection, float *current_projection, float *invViewMatrix, u_int_2 detectorPixels, float_3 sourcePosition, u_int_3 volumeVoxels, float_3 volumeSize, float t_step, int interpolation); 


