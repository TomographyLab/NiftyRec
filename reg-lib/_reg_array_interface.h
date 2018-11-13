

#include "_reg.h"

int reg_array_compute_control_point_size(int size[], int imageSize[], float gridSpacing[]);
int reg_array_bspline_initialiseControlPointGridWithAffine(float *controlPoints, float *affineTransformation, int imageSize[], float gridSpacing[]);
int reg_array_gradient_NMI_nodes(float *target, float *source, float *gradient, int *image_size, float *control_points, float gridSpacing[], int binning, int enable_gpu);
int reg_array_gradient_voxel_to_nodes(float *cp_gradient_ptr, float *gradient_ptr, float *cp_ptr, int image_size[], int cp_size[], float grid_spacing[], int enable_gpu);
int reg_array_resample_spline(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_points_size[], float *spacing, int enable_gpu);
int reg_array_image_gradient(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_points_size[], float *spacing, int enable_gpu);
int reg_array_ssd_gradient(float *ssd_gradient_ptr, float *target_ptr, float *source_ptr, float *gradient_ptr, int image_size[], int smoothing_radius[], int enable_gpu);
int reg_array_gaussian_smooth(float *image_ptr, int image_size[], float smoothing_sigma, int enable_gpu); 
int reg_array_scale_amplitude(float *image_ptr, int image_size[], float min_value, float max_value, int enable_gpu); 
int reg_array_gradient_jacobian_determinant(float *nodes_gradient_ptr, float *control_points_ptr, int image_size[], int cp_size[], float cp_spacing[], float weight, int GPU);
int reg_array_gradient_bending_energy(float *nodes_gradient_ptr, float *control_points_ptr, int image_size[], int cp_size[], float cp_spacing[], float weight, int GPU);




/* Status flags and test function for SimpleWrap Python wrapper */
extern "C" int status_success(int *status);
extern "C" int status_io_error(int *status);
extern "C" int status_initialisation_error(int *status);
extern "C" int status_parameter_error(int *status);
extern "C" int status_unhandled_error(int *status);
extern "C" int echo(int *in, int *out); 

extern "C" unsigned int REG_array_resample_image_rigid(float *inimage, float *outimage, unsigned int *image_size_x, unsigned int *image_size_y, unsigned int *image_size_z, float *translation, float *rotation, float *rotation_center, float *sform, unsigned int *enable_gpu); 
extern "C" unsigned int REG_array_d_intensity_d_space_rigid(float *inimage, float *outimage, unsigned int *image_size_x, unsigned int *image_size_y, unsigned int *image_size_z, float *translation, float *rotation, float *rotation_center, float *sform, unsigned int *enable_gpu);
extern "C" unsigned int REG_array_d_intensity_d_transformation_rigid(float *inimage, float *outimage, unsigned int *image_size_x, unsigned int *image_size_y, unsigned int *image_size_z, float *translation, float *rotation, float *rotation_center, float *sform, unsigned int *enable_gpu);

extern "C" unsigned int REG_array_d_ssd_d_transformation_rigid(float *ssd_gradient, float *source, float *target, float *translation, float *rotation, float *sform_source, float *sform_target, unsigned int *source_size, unsigned int *target_size, unsigned int *smoothing_radius, unsigned int *enable_gpu); 
extern "C" unsigned int REG_array_gaussian_smoothing(float *image, unsigned int *image_size, float *smoothing_sigma, unsigned int *enable_gpu); 