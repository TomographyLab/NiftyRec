/*
 *  _tt_array_interface.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_tt.h"

extern "C" int tt_array_project_perspective(float *attenuation, int *attenuation_size, float *projection, int *projection_size, float *image_origin, float *detector_origin, float *detector_shape, float *psf, int *psf_size, float background, int GPU);

extern "C" int tt_array_project_ray(VolumeType h_volume[], u_int_3 volume_voxels, float out_projections[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_position[], float_3 volume_size, float t_step, int use_gpu);

extern "C" int tt_array_backproject_ray(float h_projections[], u_int_2 detector_pixels, u_int n_projections, float out_backprojection[], float_2 detector_size[], float_3 detector_transl[], float_3 detector_rotat[], float_3 source_pos[], u_int_3 volume_voxels, float_3 volume_size, float t_step, int interpolation, int use_gpu);
