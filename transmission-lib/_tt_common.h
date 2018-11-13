/*
 *  _tt_common.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#ifndef _TT_COMMON_H
#define _TT_COMMON_H

#include <math.h>



#define INTERP_NEAREST 0
#define INTERP_TRILINEAR 1

typedef unsigned int  u_int;
typedef unsigned char u_char;
//typedef unsigned short VolumeType;
typedef float VolumeType;
//typedef char VolumeType;

struct float_2
{
  float x, y;
};
typedef struct float_2 float_2;

struct float_3
{
  float x, y, z;
};
typedef struct float_3 float_3;

struct float_4
{
  float x, y, z, w;
};
typedef struct float_4 float_4;

struct int_3
{
  int x, y, z;
};
typedef struct int_3 int_3;

struct u_int_3
{
  u_int x, y, z;
};
typedef struct u_int_3 u_int_3;

struct u_int_2
{
  u_int w, h;
};
typedef struct u_int_2 u_int_2;



struct mat_44{                   /** 4x4 matrix struct **/
  float m[4][4] ;
};
typedef mat_44 mat_44;

struct mat_33{                   /** 3x3 matrix struct **/
  float m[3][3] ;
};
typedef mat_33 mat_33;

extern "C" mat_44 reg_mat_44_mul(mat_44 *A, mat_44 *B);
extern "C" int create_rotation_matrix44(mat_44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z);

#endif
