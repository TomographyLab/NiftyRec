/*
 *  _et.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Harvard University, Martinos Center for Biomedical Imaging
 *  Jan. 2014.
 */


/**
* Functions for Emission Tomography
*/


#define SIZE_OF_INFO 5
#define MAX_DEVICES 16
#define ET_BLOCK_SIZE 8
#define ET_ERROR_BADGRID 2
#define eps 0.000000000001f


#include <stdbool.h>
#include "_niftyrec_memory.h"
#include "_et_common.h"
#include "_reg_tools.h"
#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_et_line_integral.h"
#include "_et_line_integral_attenuated.h"
#include "_et_accumulate.h"
#include "_et_line_backproject.h"
#include "_et_line_backproject_attenuated.h"
#include "_et_clear_accumulator.h"
#include "_et_convolve2D.h"
#include "_et_convolveSeparable2D.h"
#include "_et_convolveFFT2D.h"


#ifdef _USE_CUDA
#include "_reg_cudaCommon.h"
#include "_reg_resampling_gpu.h"
#include "_reg_affineTransformation_gpu.h"
#include "_et_line_integral_gpu.h"
#include "_et_line_integral_attenuated_gpu.h"
#include "_et_accumulate_gpu.h"
#include "_et_line_backproject_gpu.h"
#include "_et_line_backproject_attenuated_gpu.h"
#include "_et_attenuation_gradient_gpu.h"
#include "_et_clear_accumulator_gpu.h"
#include "_et_convolveFFT2D_gpu.h"
#include "_et_convolveSeparable2D_gpu.h"
#include "_pet_line_integral_compressed_gpu.h"
#include "_pet_line_backproject_compressed_gpu.h"

int et_rotate_gpu(nifti_image *sourceImage, nifti_image *resultImage, float alpha, float beta, float gamma, float center_x, float center_y, float center_z, float background, int axis_order);
int et_project_gpu(nifti_image *activity, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float backgroundAttenuation, int truncate_negative_values);
int et_backproject_gpu(nifti_image *sinogram, nifti_image *accumulator, nifti_image *psf, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values);
int et_joint_histogram_gpu(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B);
int et_project_backproject_gpu(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_alpha, float *cameras_beta, float *cameras_gamma, int truncate_negative_values);
int et_affine_gpu(nifti_image *sourceImage, nifti_image *resultImage, mat44 *affineTransformation, float background);
int et_convolve_gpu(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage);
int et_fisher_grid_gpu(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *priorfisherImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation); 
int et_gradient_attenuation_gpu(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values); 
int et_list_gpus(int *device_count_out, int *devices);
int et_set_gpu(int id);
int et_project_partial_gpu(nifti_image *activity, nifti_image *sinoImage, nifti_image *partialsumImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, nifti_image *backgroundImage, float backgroundAttenuation, int truncate_negative_values, int do_rotate_partial);
int et_reset_gpu(); 
#endif

unsigned short int et_is_block_multiple(unsigned short int size);
unsigned short int et_get_block_size(void);
int et_rotate(nifti_image *sourceImage, nifti_image *resultImage, float alpha, float beta, float gamma, float center_x, float center_y, float center_z, float background, int axis_order);
int et_project(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values);
int et_backproject(nifti_image *sinogram, nifti_image *accumulator, nifti_image *psf, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values);
int et_joint_histogram(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B);
int et_project_backproject(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_alpha, float *cameras_beta, float *cameras_gamma);
int et_affine(nifti_image *sourceImage, nifti_image *transformedImage, mat44 *affineTransformation, float background);
int et_convolve(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage);
int et_fisher_grid(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *priorfisherImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation); 
int et_gradient_attenuation(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values); 
int et_project_partial(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *partialsumImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, nifti_image *backgroundImage, float background_attenuation, int truncate_negative_values, int do_rotate_partial);


/* PET */ 
#ifdef _USE_CUDA
unsigned int PET_project_compressed_gpu(nifti_image *activityImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
                float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
                unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int truncate_negative_values, 
                unsigned int direction, unsigned int block_size, float *time_profiling); 
#endif

unsigned int PET_project_compressed_cpu(nifti_image *activityImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
                float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
                unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int truncate_negative_values,
                unsigned int direction);
#ifdef _USE_CUDA
unsigned int PET_project_compressed_gpu_test(nifti_image *activityImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal); 
#endif

#ifdef _USE_CUDA
unsigned int PET_backproject_compressed_gpu(nifti_image *backprojectionImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
                float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
                unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int direction, unsigned int block_size, float *time_profiling); 
#endif

unsigned int PET_backproject_compressed_cpu(nifti_image *backprojectionImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
                float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
                unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int direction); 

unsigned int et_spherical_phantom(nifti_image *phantomImage, float centerx, float centery, float centerz, float radius, float inner_value, float outer_value); 
unsigned int et_cylindrical_phantom(nifti_image *phantomImage, float centerx, float centery, float centerz, float radius, float length, unsigned int axis, float inner_value, float outer_value); 
unsigned int et_spheres_ring_phantom(nifti_image *phantomImage, float centerx, float centery, float centerz, float ring_radius, float min_sphere_radius, float max_sphere_radius, unsigned int N_spheres, float inner_value, float outer_value, float taper, unsigned int ring_axis); 

unsigned int tr_resample_grid_cpu(nifti_image *resampled, nifti_image *image, nifti_image *grid, float background, unsigned int interpolation_mode); 
#ifdef _USE_CUDA
unsigned int tr_resample_grid_gpu(nifti_image *resampled, nifti_image *image, nifti_image *grid, float background, unsigned int interpolation_mode); 
#endif
unsigned int ET_box_to_grid(nifti_image *image, mat44 *affine); 

#ifdef _USE_CUDA
unsigned int tr_transform_grid_gpu(nifti_image *transformed_grid, nifti_image *grid, mat44 *affine); 
#endif
unsigned int tr_transform_grid_cpu(nifti_image *transformed_grid, nifti_image *grid, mat44 *affine); 


// Compress and uncompress projection data: 
unsigned int pet_initialize_compression_structure(unsigned int *N_axial, unsigned int *N_azimuthal, unsigned int *N_u, unsigned int *N_v, int *offsets_matrix, unsigned short *locations );

unsigned int pet_compress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection); 

unsigned int pet_uncompress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection);


unsigned int pet_compress_projection_array(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, float *threshold, unsigned int *N_active, bool *active, float *data, int *offsets, unsigned short *locations, float *compressed_data); 

unsigned int pet_get_subset_sparsity(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, unsigned int *N_locations, int *offsets, unsigned short *locations, unsigned int *subsets_matrix, int *offsets_sub, unsigned short *locations_sub); 

unsigned int pet_get_subset_projection_array(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, unsigned int *N_locations, float *data, int *offsets, unsigned short *locations, int *offsets_sub, unsigned short *locations_sub, unsigned int *subsets_matrix, float *data_sub); 

