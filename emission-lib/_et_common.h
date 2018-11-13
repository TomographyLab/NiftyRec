/*
 *  _et_common.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _ET_COMMON_H
#define _ET_COMMON_H


#include<stdio.h>
#include<stdarg.h>

#include "_reg_tools.h"
#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"

//#define _VERBOSE

extern "C" int fprintf_verbose (const char *__restrict __format, ...);

#define XYZ_ROTATION 0
#define XZY_ROTATION 1
#define YXZ_ROTATION 2
#define YZX_ROTATION 3
#define ZXY_ROTATION 4
#define ZYX_ROTATION 5

int et_create_rotation_matrix(mat44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, int axis_order);
int et_create_scale_matrix(mat44 *transformationMatrix, float scale_x, float scale_y, float scale_z); 
int et_create_translation_matrix(mat44 *translationMatrix, float t_x, float t_y, float t_z); 
int et_create_affine_matrix_from_buffer(mat44 *matrix, float *buffer); 

static char _niftyrec_msg[2048]; 

#define niftyrec_success              0
#define niftyrec_error_unspecified    1
#define niftyrec_error_parameters     2
#define niftyrec_error_kernel         3
#define niftyrec_error_allocgpu       4
#define niftyrec_error_alloccpu       5
#define niftyrec_error_transfergpu    6
#define niftyrec_error_nogpubuilt     7 
#define niftyrec_error_nogpuinstalled 8 
#define STATUS_SUCCESS                0
#define STATUS_IO_ERROR               101
#define STATUS_INITIALISATION_ERROR   103
#define STATUS_PARAMETER_ERROR        104
#define STATUS_UNHANDLED_ERROR        105 
#define STATUS_MEMORY_ERROR           106

extern "C" char *niftyrec_error_msg(int status);


#endif
