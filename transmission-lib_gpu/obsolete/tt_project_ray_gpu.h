

typedef unsigned int  u_int;
typedef unsigned char u_char;
typedef unsigned short VolumeType;

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

extern "C" int tt_project_ray_array(VolumeType h_volume[], u_int_3 volume_voxels, float out_projections[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_pos[], float_3 volume_size, float t_step);



