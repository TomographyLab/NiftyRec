
#include "_reg_array_interface.h"
#include "_reg_affineTransformation.h"
#include "_reg_bspline.h"



//***************************************************************************************
/* SimpleWrap status and test functions */

extern "C" int status_success(int *status)
{
    *status = STATUS_SUCCESS; 
    return STATUS_SUCCESS; 
}
extern "C" int status_io_error(int *status)
{
    *status = STATUS_IO_ERROR; 
    return STATUS_SUCCESS; 
}
extern "C" int status_initialisation_error(int *status)
{
    *status = STATUS_INITIALISATION_ERROR; 
    return STATUS_SUCCESS; 
}
extern "C" int status_parameter_error(int *status)
{
    *status = STATUS_PARAMETER_ERROR; 
    return STATUS_SUCCESS; 
}
extern "C" int status_unhandled_error(int *status)
{
    *status = STATUS_UNHANDLED_ERROR; 
    return STATUS_SUCCESS; 
}

/* test SimpleWrap interface */
extern "C" int echo(int *in, int *out)
{
    *out = *in; 
    return STATUS_SUCCESS; 
}


//***************************************************************************************





int reg_array_compute_control_point_size(int size[], int imageSize[], float gridSpacing[])
{
	size[0]=floor(imageSize[0]/gridSpacing[0])+4;
	size[1]=floor(imageSize[1]/gridSpacing[1])+4;
	size[2]=floor(imageSize[2]/gridSpacing[2])+4;
	return 0;
}


int reg_array_bspline_initialiseControlPointGridWithAffine(float *controlPoints, float *affineTransformation, int imageSize[], float gridSpacing[])
{
	mat44 *affineTransformation_mat44 = (mat44 *)calloc(1,sizeof(mat44));
	// Initialize affine transform matrix
	
	int k=0;
	for (unsigned int i=0; i<4; i++) {
		for (unsigned int j=0; j<4; j++) {				
			affineTransformation_mat44->m[i][j] = affineTransformation[k];
			k++;
		}
	}

//fprintf(stderr,"Transformation matrix1: %f %f %f %f\n",affineTransformation_mat44->m[0][0],affineTransformation_mat44->m[0][1],affineTransformation_mat44->m[0][2],affineTransformation_mat44->m[0][3]);
//fprintf(stderr,"Transformation matrix2: %f %f %f %f\n",affineTransformation_mat44->m[1][0],affineTransformation_mat44->m[1][1],affineTransformation_mat44->m[1][2],affineTransformation_mat44->m[1][3]);
//fprintf(stderr,"Transformation matrix3: %f %f %f %f\n",affineTransformation_mat44->m[2][0],affineTransformation_mat44->m[2][1],affineTransformation_mat44->m[2][2],affineTransformation_mat44->m[2][3]);
//fprintf(stderr,"Transformation matrix4: %f %f %f %f\n",affineTransformation_mat44->m[3][0],affineTransformation_mat44->m[3][1],affineTransformation_mat44->m[3][2],affineTransformation_mat44->m[3][3]);
//fprintf(stderr,"Image size:             %d %d %d   \n",imageSize[0],imageSize[1],imageSize[2]);
//fprintf(stderr,"Grid spacing:           %f %f %f   \n",gridSpacing[0],gridSpacing[1],gridSpacing[2]);

	int control_point_size[3];
	reg_array_compute_control_point_size(control_point_size, imageSize, gridSpacing);
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, gridSpacing);
	controlPointImage->data = (float*) controlPoints; //FIXME: verify that nifti_make_new_nim() does not alloc ->data

//fprintf(stderr,"sform_code:    %d \n",controlPointImage->sform_code);

//fprintf(stderr,"Transformation matrix sto1: %f %f %f %f\n",controlPointImage->sto_xyz.m[0][0],controlPointImage->sto_xyz.m[0][1],controlPointImage->sto_xyz.m[0][2],controlPointImage->sto_xyz.m[0][3]);
//fprintf(stderr,"Transformation matrix sto2: %f %f %f %f\n",controlPointImage->sto_xyz.m[1][0],controlPointImage->sto_xyz.m[1][1],controlPointImage->sto_xyz.m[1][2],controlPointImage->sto_xyz.m[1][3]);
//fprintf(stderr,"Transformation matrix sto3: %f %f %f %f\n",controlPointImage->sto_xyz.m[2][0],controlPointImage->sto_xyz.m[2][1],controlPointImage->sto_xyz.m[2][2],controlPointImage->sto_xyz.m[2][3]);
//fprintf(stderr,"Transformation matrix sto4: %f %f %f %f\n",controlPointImage->sto_xyz.m[3][0],controlPointImage->sto_xyz.m[3][1],controlPointImage->sto_xyz.m[3][2],controlPointImage->sto_xyz.m[3][3]);

//fprintf(stderr,"Transformation matrix qto1: %f %f %f %f\n",controlPointImage->qto_xyz.m[0][0],controlPointImage->qto_xyz.m[0][1],controlPointImage->qto_xyz.m[0][2],controlPointImage->qto_xyz.m[0][3]);
//fprintf(stderr,"Transformation matrix qto2: %f %f %f %f\n",controlPointImage->qto_xyz.m[1][0],controlPointImage->qto_xyz.m[1][1],controlPointImage->qto_xyz.m[1][2],controlPointImage->qto_xyz.m[1][3]);
//fprintf(stderr,"Transformation matrix qto3: %f %f %f %f\n",controlPointImage->qto_xyz.m[2][0],controlPointImage->qto_xyz.m[2][1],controlPointImage->qto_xyz.m[2][2],controlPointImage->qto_xyz.m[2][3]);
//fprintf(stderr,"Transformation matrix qto4: %f %f %f %f\n",controlPointImage->qto_xyz.m[3][0],controlPointImage->qto_xyz.m[3][1],controlPointImage->qto_xyz.m[3][2],controlPointImage->qto_xyz.m[3][3]);

	if(reg_bspline_initialiseControlPointGridWithAffine(affineTransformation_mat44, controlPointImage)) return 1;

	//Free
	free_nifti_image_except_data(controlPointImage);

	return 0;
}



int reg_array_gradient_NMI_nodes(float *target, float *source, float *gradient, int *image_size, float *control_points, float gridSpacing[], int binning, int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(source);
	
	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);
	targetImage->data = (float *)(target);

	//Create control points nifti image
	int control_point_size[3];
	reg_array_compute_control_point_size(control_point_size, image_size, gridSpacing);
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, gridSpacing);
	controlPointImage->data = (float*) (control_points); //FIXME: verify that nifti_make_new_nim() does not alloc ->data

	//Create control points gradient nifti image
	nifti_image *gradientImage=NULL;
	gradientImage = nifti_copy_nim_info(controlPointImage);
	gradientImage->datatype = NIFTI_TYPE_FLOAT32;
	gradientImage->nbyper = sizeof(float);
	gradientImage->data = (float*) (gradient);

	//Compute NMI gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gradient_NMI_nodes_gpu(targetImage, sourceImage, controlPointImage, gradientImage, binning);
	else
		status = reg_gradient_NMI_nodes(targetImage, sourceImage, controlPointImage, gradientImage, binning);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_NMI_nodes_array: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	status = reg_gradient_NMI_nodes(targetImage, sourceImage, controlPointImage, gradientImage, binning);
	#endif

	//Free
	free_nifti_image_except_data(targetImage);
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(controlPointImage);
	free_nifti_image_except_data(gradientImage);

	return status;
}


int reg_array_gradient_voxel_to_nodes(float *control_points_gradient, float *gradient, float *control_points, int image_size[], int cp_size[], float grid_spacing[], int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *gradientImage = reg_initialize_deformation_field(image_size, image_spacing);
	gradientImage->data = (float *)(gradient);

	//Create control points nifti image
	int control_point_size[3];
	reg_array_compute_control_point_size(control_point_size, image_size, grid_spacing);
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, grid_spacing);
	controlPointImage->data = (float*) (control_points); //FIXME: verify that nifti_make_new_nim() does not alloc ->data

fprintf(stderr,"Size:     %d %d %d %d %d \n",gradientImage->nx,gradientImage->ny,gradientImage->nz,gradientImage->nt,gradientImage->nu); 
fprintf(stderr,"Delta:    %f %f %f %f %f \n",gradientImage->dx,gradientImage->dy,gradientImage->dz,gradientImage->dt,gradientImage->du); 
fprintf(stderr,"Data:     %f %f %f %f \n",((float*)gradientImage->data)[0],((float*)gradientImage->data)[1],((float*)gradientImage->data)[2],((float*)gradientImage->data)[3]); 

	//Create control points gradient nifti image
	nifti_image *cpGradientImage=NULL;
	cpGradientImage = nifti_copy_nim_info(controlPointImage);
	cpGradientImage->datatype = NIFTI_TYPE_FLOAT32;
	cpGradientImage->nbyper = sizeof(float);
	cpGradientImage->data = (float*) control_points_gradient;

	//Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gradient_voxel_to_nodes_gpu(cpGradientImage, gradientImage, controlPointImage);
	else
		status = reg_gradient_voxel_to_nodes(cpGradientImage, gradientImage, controlPointImage);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_NMI_nodes_array: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	status = reg_gradient_voxel_to_nodes(cpGradientImage, gradientImage, controlPointImage);
	#endif

	//Free
	free_nifti_image_except_data(gradientImage);
	free_nifti_image_except_data(cpGradientImage);
	free_nifti_image_except_data(controlPointImage);

	return status;
}



int reg_array_resample_spline(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_point_size[], float *spacing, int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(image_ptr);

	nifti_image *resultImage = reg_initialize_image(image_size, image_spacing);
	resultImage->data = (float *)(outimage_ptr);

	//Create control points nifti image
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, spacing);
	controlPointImage->data = (float*) nodes_ptr; //FIXME: verify that nifti_make_new_nim() does not alloc ->data

	//Resample image
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_resample_spline_gpu(resultImage, sourceImage, controlPointImage);
	else
		status = reg_resample_spline(resultImage, sourceImage, controlPointImage);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_NMI_nodes_array: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	status = reg_resample_spline(resultImage, sourceImage, controlPointImage);
	#endif	

	//Free 
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(resultImage);
	free_nifti_image_except_data(controlPointImage);

	return status;
}



int reg_array_image_gradient(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_points_size[], float *spacing, int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(image_ptr);

	nifti_image *resultImage = reg_initialize_deformation_field(image_size, image_spacing);
	resultImage->data = (float *)(outimage_ptr);

	//Create control points nifti image
	nifti_image *controlPointImage = reg_initialize_control_points(control_points_size, spacing);
	controlPointImage->data = (float*) nodes_ptr; 

	//Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_image_gradient_gpu(resultImage, sourceImage, controlPointImage);
	else
		status = reg_image_gradient(resultImage, sourceImage, controlPointImage);
	#else
	if (GPU)
		fprintf(stderr, "reg_image_gradient: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = reg_image_gradient(resultImage, sourceImage, controlPointImage);
	#endif	
	
	//Free 
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(resultImage);
	free_nifti_image_except_data(controlPointImage);

	return status;
}


int reg_array_ssd_gradient(float *ssd_gradient_ptr, float *target_ptr, float *source_ptr, float *gradient_ptr, int image_size[], int smoothingRadius[], int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);
	targetImage->data = (float *)(target_ptr);

	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(source_ptr);

	nifti_image *gradientImage = reg_initialize_deformation_field(image_size, image_spacing);
	gradientImage->data = (float *)(gradient_ptr);

	nifti_image *ssdGradientImage = reg_initialize_deformation_field(image_size, image_spacing);
	ssdGradientImage->data = (float *)(ssd_gradient_ptr);

        // Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_ssd_gradient_gpu(ssdGradientImage, targetImage, sourceImage, gradientImage, smoothingRadius);
	else
		status = reg_ssd_gradient(ssdGradientImage, targetImage, sourceImage, gradientImage, smoothingRadius);
	#else
	if (GPU)
		fprintf(stderr, "reg_ssd_gradient: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = reg_ssd_gradient(ssdGradientImage, targetImage, sourceImage, gradientImage, smoothingRadius);
	#endif	

	//Free 
	free_nifti_image_except_data(targetImage);
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(gradientImage);
	free_nifti_image_except_data(ssdGradientImage);

        return status;
}


int reg_array_gaussian_smooth(float *image_ptr, int image_size[], float smoothing_sigma, int GPU)
{
	int status = 1;
	// Allocate image
	float image_spacing[3] = {1,1,1};
	nifti_image *niftiImage = reg_initialize_image(image_size, image_spacing);
	niftiImage->data = (float *)(image_ptr);

        // Smooth
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gaussian_smooth_gpu(niftiImage, smoothing_sigma);
	else
		status = reg_gaussian_smooth(niftiImage, smoothing_sigma);
	#else
	if (GPU)
		fprintf(stderr, "reg_gaussian_smoothing: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = reg_gaussian_smooth(niftiImage, smoothing_sigma);
	#endif	

	// Free
	free_nifti_image_except_data(niftiImage);

	return status;
}


int reg_array_scale_amplitude(float *image_ptr, int image_size[], float min_value, float max_value, int GPU)
{
	int status = 1;
	// Allocate image 
	float image_spacing[3] = {1,1,1};
	nifti_image *niftiImage = reg_initialize_image(image_size, image_spacing);
	niftiImage->data = (float *)(image_ptr);

        // Smooth
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_scale_amplitude_gpu(niftiImage, min_value, max_value);
	else
		status = reg_scale_amplitude(niftiImage, min_value, max_value);
	#else
	if (GPU)
		fprintf(stderr, "reg_scale_amplitude: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = reg_scale_amplitude(niftiImage, min_value, max_value);
	#endif	

	// Free
	free_nifti_image_except_data(niftiImage);
	
	return status;
}


int reg_array_gradient_jacobian_determinant(float *nodes_gradient_ptr, float *control_points_ptr, int image_size[], int cp_size[], float cp_spacing[], float weight, int GPU)
{
	int status = 1;
	float image_spacing[3] = {1,1,1};

	//Create control points nifti image
	nifti_image *controlPointImage = reg_initialize_control_points(cp_size, cp_spacing);
	controlPointImage->data = (float*) control_points_ptr; 

	//Create nifti image for gradient (result)
	nifti_image *nodesGradientImage = reg_initialize_control_points(cp_size, cp_spacing);
	nodesGradientImage->data = (float*) nodes_gradient_ptr; 

	//Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gradient_jacobian_determinant_gpu(nodesGradientImage, controlPointImage, image_size, image_spacing, weight);
	else
		status = reg_gradient_jacobian_determinant(nodesGradientImage, controlPointImage, image_size, image_spacing, weight);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_jacobian_determinant: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = reg_gradient_jacobian_determinant(nodesGradientImage, controlPointImage, image_size, image_spacing, weight);
	#endif	

	//Free 
	free_nifti_image_except_data(controlPointImage);
	free_nifti_image_except_data(nodesGradientImage);

	return status;
}


int reg_array_gradient_bending_energy(float *nodes_gradient_ptr, float *control_points_ptr, int image_size[], int cp_size[], float cp_spacing[], float weight, int GPU)
{
	int status = 1;
	float image_spacing[3] = {1,1,1};

	//Create control points nifti image
	nifti_image *controlPointImage = reg_initialize_control_points(cp_size, cp_spacing);
	controlPointImage->data = (float*) control_points_ptr; 

	//Create nifti image for gradient (result)
	nifti_image *nodesGradientImage = reg_initialize_control_points(cp_size, cp_spacing);
	nodesGradientImage->data = (float*) nodes_gradient_ptr; 

	//Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gradient_bending_energy_gpu(nodesGradientImage, controlPointImage, image_size, image_spacing, weight);
	else
		status = reg_gradient_bending_energy(nodesGradientImage, controlPointImage, image_size, image_spacing, weight);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_bending_energy: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = reg_gradient_bending_energy(nodesGradientImage, controlPointImage, image_size, image_spacing, weight);
	#endif	

	//Free 
	free_nifti_image_except_data(controlPointImage);
	free_nifti_image_except_data(nodesGradientImage);

	return status;
}











extern "C" unsigned int REG_array_resample_image_rigid(float *inimage, float *outimage, unsigned int *image_size_x, unsigned int *image_size_y, unsigned int *image_size_z, float *translation, float *rotation, float *rotation_center, float *sform, unsigned int *GPU) 
{
    int status=STATUS_SUCCESS; 
    
    // Define nifti structure for the input image 
    int dim[8]; 
    dim[0]    = 3; 
    dim[1]    = *image_size_x; 
    dim[2]    = *image_size_y; 
    dim[3]    = *image_size_z; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *inputImage = NULL; 
    inputImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    inputImage->data = (float *)(inimage); 

    inputImage->pixdim[0]  = 1;    
    inputImage->qform_code = 0; 
    inputImage->sform_code = 1; 
    inputImage->sto_xyz.m[0][0]=sform[0];     inputImage->sto_xyz.m[0][1]=sform[1];    inputImage->sto_xyz.m[0][2]=sform[2];     inputImage->sto_xyz.m[0][3]=sform[3];
    inputImage->sto_xyz.m[1][0]=sform[4];     inputImage->sto_xyz.m[1][1]=sform[5];    inputImage->sto_xyz.m[1][2]=sform[6];     inputImage->sto_xyz.m[1][3]=sform[7];
    inputImage->sto_xyz.m[2][0]=sform[8];     inputImage->sto_xyz.m[2][1]=sform[9];    inputImage->sto_xyz.m[2][2]=sform[10];    inputImage->sto_xyz.m[2][3]=sform[11]; 
    inputImage->sto_xyz.m[3][0]=sform[12];    inputImage->sto_xyz.m[3][1]=sform[13];   inputImage->sto_xyz.m[3][2]=sform[14];    inputImage->sto_xyz.m[3][3]=sform[15];
    inputImage->sto_ijk = nifti_mat44_inverse(inputImage->sto_xyz); 

    // Define nifti structure for the resampled image 
    dim[0]    = 3; 
    dim[1]    = *image_size_x; 
    dim[2]    = *image_size_y; 
    dim[3]    = *image_size_z; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *resampledImage = NULL; 
    resampledImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    resampledImage->data = (float *)(outimage);
 
    resampledImage->pixdim[0]  = 1;    
    resampledImage->qform_code = 0; 
    resampledImage->sform_code = 1; 
    resampledImage->sto_xyz.m[0][0]=sform[0];     resampledImage->sto_xyz.m[0][1]=sform[1];    resampledImage->sto_xyz.m[0][2]=sform[2];     resampledImage->sto_xyz.m[0][3]=sform[3];
    resampledImage->sto_xyz.m[1][0]=sform[4];     resampledImage->sto_xyz.m[1][1]=sform[5];    resampledImage->sto_xyz.m[1][2]=sform[6];     resampledImage->sto_xyz.m[1][3]=sform[7];
    resampledImage->sto_xyz.m[2][0]=sform[8];     resampledImage->sto_xyz.m[2][1]=sform[9];    resampledImage->sto_xyz.m[2][2]=sform[10];    resampledImage->sto_xyz.m[2][3]=sform[11]; 
    resampledImage->sto_xyz.m[3][0]=sform[12];    resampledImage->sto_xyz.m[3][1]=sform[13];   resampledImage->sto_xyz.m[3][2]=sform[14];    resampledImage->sto_xyz.m[3][3]=sform[15];
    resampledImage->sto_ijk = nifti_mat44_inverse(resampledImage->sto_xyz); 

    // Define the affine transformation matrix for the rigid transformation 
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *rotation_m    = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *m             = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_translation_matrix(translation_m, translation[0], translation[1], translation[2]); 
    et_create_rotation_matrix(rotation_m, rotation[0], rotation[1], rotation[2], rotation_center[0], rotation_center[1], rotation_center[2] , XYZ_ROTATION); 
	*m = reg_mat44_mul(rotation_m, translation_m);

    // Resample the image 
	#ifdef _USE_CUDA
	if (*GPU)
		status = REG_resample_image_rigid_gpu(resampledImage, inputImage, m); 
	else
		status = REG_resample_image_rigid_cpu(resampledImage, inputImage, m); 
	#else
	if (*GPU) {
		fprintf(stderr, "REG_resample_image_rigid_cpu: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	}
	status = REG_resample_image_rigid_cpu(resampledImage, inputImage, m); 
	#endif	    

    // Free
    free(translation_m); 
    free(rotation_m); 
    free_nifti_image_except_data(inputImage);
	free_nifti_image_except_data(resampledImage);

    return status; 
}





extern "C" unsigned int REG_array_d_intensity_d_space_rigid(float *inimage, float *outimage, unsigned int *image_size_x, unsigned int *image_size_y, unsigned int *image_size_z, float *translation, float *rotation, float *rotation_center, float *sform, unsigned int *GPU)
{
	int status;
    int image_size[3]; 
    image_size[0]=*image_size_x;
    image_size[1]=*image_size_y;
    image_size[2]=*image_size_z;
    float image_spacing[3] = {1,1,1};
    	
	// Allocate source and target nifti image 
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(inimage);

    // Set sform  .. FIXME: do it

	nifti_image *resultImage = reg_initialize_deformation_field(image_size, image_spacing);
	resultImage->data = (float *)(outimage);

    // Create affine matrix for rigid transformation: 
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *rotation_m    = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *affine        = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_translation_matrix(translation_m, translation[0], translation[1], translation[2]); 
    et_create_rotation_matrix(rotation_m, rotation[0], rotation[1], rotation[2], rotation_center[0], rotation_center[1], rotation_center[2] , XYZ_ROTATION); 
	*affine = reg_mat44_mul(rotation_m, translation_m);

	// Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = REG_d_intensity_d_space_rigid_gpu(resultImage, sourceImage, affine);
	else
		status = REG_d_intensity_d_space_rigid_cpu(resultImage, sourceImage, affine);
	#else
	if (GPU)
		fprintf(stderr, "REG_array_d_intensity_d_space_rigid: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = REG_d_intensity_d_space_rigid_cpu(resultImage, sourceImage, affine);
	#endif	
	
	// Free 
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(resultImage);
    free(affine); 
    free(translation_m); 
    free(rotation_m); 
    
	return status;
}



extern "C" unsigned int REG_array_d_intensity_d_transformation_rigid(float *inimage, float *outimage, unsigned int *image_size_x, unsigned int *image_size_y, unsigned int *image_size_z, float *translation, float *rotation, float *rotation_center, float *sform, unsigned int *GPU)
{
	int status;
    int image_size[4]; 
    image_size[0]=*image_size_x;
    image_size[1]=*image_size_y;
    image_size[2]=*image_size_z;
    image_size[3]=6; 
    float image_spacing[3] = {1,1,1};
    	
	// Allocate source and target nifti image 
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing); 
	sourceImage->data = (float *)(inimage);

    // Set sform  .. FIXME: do it

	nifti_image *resultImage = reg_initialize_image_4D(image_size, image_spacing); 
	resultImage->data = (float *)(outimage); 

    // Create affine matrix for rigid transformation: 
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *rotation_m    = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *affine        = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_translation_matrix(translation_m, translation[0], translation[1], translation[2]); 
    et_create_rotation_matrix(rotation_m, rotation[0], rotation[1], rotation[2], rotation_center[0], rotation_center[1], rotation_center[2] , XYZ_ROTATION); 
	*affine = reg_mat44_mul(rotation_m, translation_m);

	// Compute gradients
	#ifdef _USE_CUDA
	if (GPU)
		status = REG_d_intensity_d_transformation_rigid_gpu(resultImage, sourceImage, affine);
	else
		status = REG_d_intensity_d_transformation_rigid_cpu(resultImage, sourceImage, affine);
	#else
	if (GPU)
		fprintf(stderr, "REG_array_d_intensity_d_space_rigid: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
		status = REG_d_intensity_d_transformation_rigid_cpu(resultImage, sourceImage, affine);
	#endif	
	
	// Free 
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(resultImage);
    free(affine); 
    free(translation_m); 
    free(rotation_m); 
    
	return status;
}




extern "C" unsigned int REG_array_d_ssd_d_transformation_rigid(float *ssd_gradient, float *source, float *target, float *translation, float *rotation, float *sform_source, float *sform_target, unsigned int *source_size, unsigned int *target_size, unsigned int *smoothing_radius, unsigned int *enable_gpu)
{
    int status = 0; 
    return status; 
}


extern "C" unsigned int REG_array_gaussian_smoothing(float *image, unsigned int *image_size, float *smoothing_sigma, unsigned int *enable_gpu)
{
    int status = 0; 
    return status; 
}





