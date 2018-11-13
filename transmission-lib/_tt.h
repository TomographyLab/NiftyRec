/*
 *  _tt.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_reg_tools.h"
#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_et_line_integral.h"
#include "_et_accumulate.h"
#include "_et_line_backproject.h"
#include "_et_clear_accumulator.h"
//#include "_tt_project_ray.h"
#include "_tt_common.h"
#include "_tt_line_project_ray_cpu.h"
#include "_tt_line_backproject_ray_cpu.h"

#ifdef _USE_CUDA
#include "_reg_cudaCommon.h"
#include "_reg_resampling_gpu.h"
#include "_reg_affineTransformation_gpu.h"
#include "_et_line_integral_gpu.h"
#include "_et_accumulate_gpu.h"
#include "_et_line_backproject_gpu.h"
#include "_et_clear_accumulator_gpu.h"
#include "_et_convolveFFT2D_gpu.h"
//#include "_tt_perspective_gpu.h"
#include "_tt_line_project_ray_gpu.h"
#include "_tt_line_backproject_ray_gpu.h"
#endif


int tt_backproject_ray(float h_projections[], u_int_2 detector_pixels, u_int n_projections, float out_backprojection[], float_2 detector_size[], float_3 detector_transl[], float_3 detector_rotat[], float_3 source_pos[], u_int_3 volume_voxels, float_3 volume_size, float t_step, int interpolation, int use_gpu);

int tt_project_ray(VolumeType h_volume[], u_int_3 volume_voxels, float out_projections[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_position[], float_3 volume_size, float t_step, int use_gpu);

int tt_project_perspective(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background);

#ifdef _USE_CUDA
int tt_project_perspective_gpu(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background);
#endif
