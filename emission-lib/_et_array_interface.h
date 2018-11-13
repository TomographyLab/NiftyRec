/*
 *  _et_array_interface.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Harvard University, Martinos Center for Biomedical Imaging
 *  Jan. 2014.
 */



/**
* Plain C-types interface to NiftyRec
*/

#include "_et.h"


/* Status flags and test function for SimpleWrap Python wrapper */
extern "C" int status_success(int *status);
extern "C" int status_io_error(int *status);
extern "C" int status_initialisation_error(int *status);
extern "C" int status_parameter_error(int *status);
extern "C" int status_unhandled_error(int *status);
extern "C" int echo(int *in, int *out); 


/* Emission tomography common */
extern "C" int et_array_affine(float *image_ptr, int *image_size, float *transformed_image_ptr, int *transformed_image_size, float *affine_ptr, int *affine_size, float background, int GPU);
extern "C" int et_array_rotate(float *image_ptr, int *size_ptr, float *rotated_image_ptr, float *angles_ptr, float *centers_ptr, float background, int axis_order, int GPU);
extern "C" int et_array_convolve(float *image_ptr, int *image_size, float *out_ptr, int *out_size, float *psf_ptr, int *psf_size, int GPU);
extern "C" int et_array_histogram_weighted(float *inputdata_ptr, float *weights_ptr, float *histogram_ptr, int N, int N_classes, int N_bins, double min_value, double max_value); 

extern "C" int et_array_list_gpus(int *gpu_count, int *gpus_info_array);
extern "C" int et_array_set_gpu(int id);
extern "C" int et_array_isinstalled();
extern "C" int et_array_reset_gpu();
extern "C" unsigned short int et_array_is_block_multiple(unsigned short int size); 
extern "C" unsigned short int et_array_get_block_size(void); 

/* SPECT */
extern "C" int SPECT_project_parallelholes(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float *Background, float *Background_attenuation, int *Gpu, int *Truncate_negative_values); 
extern "C" int SPECT_backproject_parallelholes(float *sino, int *sino_size, float *bkpr, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float *Background, float *Background_attenuation, int *Gpu, int *Truncate_negative_values); 

extern "C" int et_array_project_partial(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *partialsum, int *partialsum_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float *background_image, float background_attenuation, int GPU, int truncate_negative_values, int do_rotate_partial);
extern "C" int et_array_gradient_attenuation(float *sino, int *sino_size, float *activity, int *activity_size, float *gradient, int *gradient_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU, int truncate_negative_values);

extern "C" int et_array_fisher_grid(float *activity_ptr, int *activity_size, float *cameras_ptr, int *cameras_size, float *psf_ptr, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation_ptr, int *attenuation_size, float epsilon, float background, float background_attenuation, int enable_gpu);
extern "C" int et_array_fisher_grid_projection(float *sinogram_ptr, int *sinogram_size, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU);

extern "C" int et_array_calculate_size_psf(unsigned int *psf_size_x, unsigned int *psf_size_y, float fwhm_pixels_dist0, float sensitivity0, float dist0, float fwhm_pixels_dist1, float sensitivity1, float dist1); 
extern "C" int et_array_make_psf(float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float fwhm_pixels_dist0, float sensitivity0, float dist0, float fwhm_pixels_dist1, float sensitivity1, float dist1, unsigned int N_psf_planes); 
extern "C" int et_array_make_cameras(float *cameras_data, float firstcamera_deg, float lastcamera_deg, unsigned int n_cameras, unsigned int rotation_axis);
extern "C" int et_array_osem(float *activity_data, unsigned int size_x, unsigned int size_y, unsigned int subset_order, float *sinogram_data, int n_cameras, float firstcamera, float lastcamera, unsigned int rotation_axis, int iterations, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
extern "C" int et_array_osem_step(float *activity_data, int size_x, int size_y, unsigned int subset_order, float *sinogram_data, int n_cameras, float *cameras_data, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
extern "C" int et_array_mlem(float *activity_data, unsigned int size_x, unsigned int size_y, float *sinogram_data, int n_cameras, float firstcamera, float lastcamera, unsigned int rotation_axis, int iterations, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
extern "C" int et_array_mlem_step(float *activity_data, int size_x, int size_y, float *sinogram_data, int n_cameras, float *cameras_data, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);


/* Transmission Tomography */
extern "C" int tt_array_project(); 
extern "C" int tt_array_backproject(); 




/* PET */
extern "C" int PET_project(int *parameter); 
extern "C" unsigned int PET_backproject(void); 

extern "C" int PET_project_compressed(
                           float *projection, 
                           float *activity, 
                           unsigned int *N_activity_x, unsigned int *N_activity_y, unsigned int *N_activity_z, 
                           float *activity_size_x, float *activity_size_y, float *activity_size_z, 
                           float *T_activity_x, float *T_activity_y, float *T_activity_z, 
                           float *R_activity_x, float *R_activity_y, float *R_activity_z, 
                           float *attenuation, 
                           unsigned int *N_attenuation_x, unsigned int *N_attenuation_y, unsigned int *N_attenuation_z, 
                           float *attenuation_size_x, float *attenuation_size_y, float *attenuation_size_z, 
                           float *T_attenuation_x, float *T_attenuation_y, float *T_attenuation_z, 
                           float *R_attenuation_x, float *R_attenuation_y, float *R_attenuation_z, 
                           unsigned int *N_axial, unsigned int *N_azimuthal, 
                           float *angles, 
                           unsigned int *N_u, unsigned int *N_v, 
                           float *size_u, float *size_v, 
                           unsigned int *N_locations, 
                           int *offsets, 
                           unsigned short *locations, 
                           unsigned int *active, 
                           unsigned int *N_samples, float *sample_step, float *background, float *background_attenuation, unsigned int *truncate_negative_values, 
                           unsigned int *use_gpu,
                           unsigned int *direction, unsigned int *block_size, float *time_profiling); 

extern "C" int PET_project_compressed_test(
                           float *projection, 
                           float *activity, 
                           unsigned int *N_activity_x, unsigned int *N_activity_y, unsigned int *N_activity_z, 
                           float *attenuation, 
                           unsigned int *N_attenuation_x, unsigned int *N_attenuation_y, unsigned int *N_attenuation_z, 
                           unsigned int *N_axial, unsigned int *N_azimuthal, 
                           unsigned int *N_locations, 
                           int *offsets, 
                           unsigned short *locations,
                           unsigned int *active ); 

extern "C" unsigned int PET_backproject_compressed(
                           float *backprojection, 
                           unsigned int *N_activity_x, unsigned int *N_activity_y, unsigned int *N_activity_z, 
                           float *activity_size_x, float *activity_size_y, float *activity_size_z, 
                           float *T_activity_x, float *T_activity_y, float *T_activity_z, 
                           float *R_activity_x, float *R_activity_y, float *R_activity_z, 
                           float *attenuation, 
                           unsigned int *N_attenuation_x, unsigned int *N_attenuation_y, unsigned int *N_attenuation_z, 
                           float *attenuation_size_x, float *attenuation_size_y, float *attenuation_size_z, 
                           float *T_attenuation_x, float *T_attenuation_y, float *T_attenuation_z, 
                           float *R_attenuation_x, float *R_attenuation_y, float *R_attenuation_z, 
                           unsigned int *N_axial, unsigned int *N_azimuthal, 
                           float *angles,  
                           unsigned int *N_u, unsigned int *N_v, 
                           float *size_u, float *size_v, 
                           unsigned int *N_locations, 
                           int *offsets, 
                           unsigned short *locations, 
                           unsigned int *active, 
                           float *projection_data, 
                           unsigned int *use_gpu,
                           unsigned int *N_samples,
                           float *sample_step,
                           float *background_activity,
                           float *background_attenuation,
                           unsigned int *direction,
                           unsigned int *block_size, float *time_profiling ); 


extern "C" unsigned int ET_spherical_phantom(float *image_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *sizex, float *sizey, float *sizez, float *centerx, float *centery, float *centerz, float *radius, float *inner_value, float *outer_value); 
extern "C" unsigned int ET_cylindrical_phantom(float *image_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *sizex, float *sizey, float *sizez, float *centerx, float *centery, float *centerz, float *radius, float *length, unsigned int *axis, float *inner_value, float *outer_value); 
extern "C" unsigned int ET_spheres_ring_phantom(float *image_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *sizex, float *sizey, float *sizez, float *centerx, float *centery, float *centerz, float *ring_radius, float *min_sphere_radius, float *max_sphere_radius, unsigned int *N_spheres, float *inner_value, float *outer_value, float *taper, unsigned int *ring_axis);

extern "C" unsigned int TR_resample_grid(float *resampled_image_array, float *image_array, float *affine, float *grid_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode); 
extern "C" unsigned int TR_resample_box(float *resampled_image_array, float *image_array, float *affine, float *min_coords, float *max_coords, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode); 
extern "C" unsigned int TR_gradient_grid(float *gradient_array, float *image_array, float *affine, float *grid, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode);
extern "C" unsigned int TR_gradient_box(float *gradient_array, float *image_array, float *affine, float *min_coords, float *max_coords, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode); 

extern "C" unsigned int TR_grid_from_box_and_affine(float *grid, float *affine_box2grid, float *min_x, float *min_y, float *min_z, float *max_x, float *max_y, float *max_z, unsigned int *n_x, unsigned int * n_y, unsigned int *n_z); 

extern "C" unsigned int TR_transform_grid(float *transformed_grid_array, float *grid_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *affine_from_grid, unsigned int *use_gpu); 

// Compress and uncompress projection data 
extern "C" unsigned int PET_initialize_compression_structure(unsigned int *N_axial, unsigned int *N_azimuthal, unsigned int *N_u, unsigned int *N_v, int *offsets_matrix, unsigned short *locations ); 
extern "C" unsigned int PET_compress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection); 
extern "C" unsigned int PET_uncompress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection); 

extern "C" unsigned int PET_compress_projection_array(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, float *threshold, unsigned int *N_active, bool *active, float *data, int *offsets, unsigned short *locations, float *compressed_data); 
extern "C" unsigned int PET_get_subset_sparsity(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, unsigned int *N_locations, int *offsets, unsigned short *locations, unsigned int *subsets_matrix, int *offsets_sub, unsigned short *locations_sub); 
extern "C" unsigned int PET_get_subset_projection_array(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, unsigned int *N_locations, float *data, int *offsets, unsigned short *locations, int *offsets_sub, unsigned short *locations_sub, unsigned int *subsets_matrix, float *data_sub);

