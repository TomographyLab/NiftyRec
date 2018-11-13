
#define PrecisionTYPE float

//extern "C" int seg_array_initialise(int numb_classes, int size_x, int size_y, int size_z, int n_dimensions, PrecisionTYPE *image_data);
extern "C" int seg_array_initialise_fromfile(int numb_classes, char* filename, char *mask_filename);

extern "C" int seg_array_set_MRF_strength(double strength);
extern "C" int seg_array_get_MRF_strength(double *MRF_strength);
extern "C" int seg_array_set_regularisation_covariance(double parameter);
extern "C" int seg_array_get_regularisation_covariance(double *parameter);
extern "C" int seg_array_set_biasfield_parameters(unsigned int order, double threshold);
extern "C" int seg_array_get_biasfield_parameters(int *enabled, unsigned int *order, double *threshold);
extern "C" int seg_array_get_biasfield_size(int* size_x, int* size_y, int* size_z);
extern "C" int seg_array_get_biasfield(PrecisionTYPE* biasfield_data);
extern "C" int seg_array_save_biasfield(char *filename);

extern "C" int seg_array_set_mask_fromfile(char *filename);
extern "C" int seg_array_set_priors_fromfiles(int n_classes, char **filenames);

extern "C" int seg_array_step(int n);
extern "C" int seg_array_step_Gaussian();
extern "C" int seg_array_step_Expectation();
extern "C" int seg_array_step_Maximization();
extern "C" int seg_array_step_MRF();
extern "C" int seg_array_step_BiasField();
extern "C" int seg_array_step_PriorWeight();

extern "C" int seg_array_get_image_size(int* size_x, int* size_y, int* size_z, int* n_volumes);
extern "C" int seg_array_get_image(PrecisionTYPE* image_data);
extern "C" int seg_array_set_input_image(PrecisionTYPE* image_data);
extern "C" int seg_array_get_segmentation_size(int* size_x, int* size_y, int* size_z, int* n_classes);
extern "C" int seg_array_get_segmentation(PrecisionTYPE* segmentation_data);
extern "C" int seg_array_set_segmentation(PrecisionTYPE* segmentation_data);
extern "C" int seg_array_save_segmentation(char *filename);
extern "C" int seg_array_get_mask(PrecisionTYPE* mask_data);
extern "C" int seg_array_save_mask(char *filename);
extern "C" int seg_array_get_priors(PrecisionTYPE* prior_data);
extern "C" int seg_array_save_priors(char *filename);
extern "C" int seg_array_get_Gaussian_parameters();
extern "C" int seg_array_get_loglikelihood(double *loglikelihood);
extern "C" int seg_array_get_mean_variance(double *mean, double *variance);
extern "C" int seg_array_set_mean(float *mean);
extern "C" int seg_array_set_variance(float *variance);

extern "C" int seg_array_erode(float *image_ptr, int *image_size, float *out_ptr, int iterations);

extern "C" int seg_array_cleanup();
