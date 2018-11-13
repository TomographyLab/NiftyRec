
#ifndef _TT_BACKPROJECT_RAY_GPU_H
#define _TT_BACKPROJECT_RAY_GPU_H


typedef unsigned int  u_int;
typedef unsigned char u_char;

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

extern "C" int tt_backproject_ray_array(float h_projections[], u_int_3 volume_voxels, float out_backprojection[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_pos[], float_3 volume_size, float t_step);

#endif
