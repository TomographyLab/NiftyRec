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

#include "_et_array_interface.h"
#include "_et_common.h"
#define PI 3.141592653589793f



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
/* Emission tomography common */

extern "C" int et_array_affine(float *image_ptr, int *image_size, float *transformed_image_ptr, int *transformed_image_size, float *affine_ptr, int *affine_size, float background, int GPU)
{
	int status;

	// Allocate source image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *sourceImage;
        sourceImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sourceImage->data = static_cast<void *>(image_ptr);
	
	// Allocate the result image
	nifti_image *transformedImage = nifti_copy_nim_info(sourceImage);
        transformedImage->data = static_cast<void *>(transformed_image_ptr);

        // Allocate and initialise Affine transformation matrix
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        //FIXME

	// Rotate
	#ifdef _USE_CUDA
	if(GPU)
	    status = et_affine_gpu(sourceImage, transformedImage, affineTransformation, background);
	else
	    status = et_affine(sourceImage, transformedImage, affineTransformation, background);
      	#else
      	    status = et_affine(sourceImage, transformedImage, affineTransformation, background);
      	#endif
      	
	//Free
	if( sourceImage->fname != NULL ) free(sourceImage->fname) ;
	if( sourceImage->iname != NULL ) free(sourceImage->iname) ;
	(void)nifti_free_extensions( sourceImage ) ;
	free(sourceImage) ;
	
	if( transformedImage->fname != NULL ) free(transformedImage->fname) ;
	if( transformedImage->iname != NULL ) free(transformedImage->iname) ;
	(void)nifti_free_extensions( transformedImage ) ;
	free(transformedImage) ;

        free(affineTransformation);

	return status;
}



extern "C" int et_array_rotate(float *image, int *size, float *rotated_image, float *angles, float *centers, float background, int axis_order, int GPU)
{
	int status;

	// Allocate source image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = size[0];
	dim[2]    = size[1];
	dim[3]    = size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *sourceImage;
        sourceImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sourceImage->data = static_cast<void *>(image);
	
        // Allocate the result image
	nifti_image *resultImage = nifti_copy_nim_info(sourceImage);
        resultImage->data = static_cast<void *>(rotated_image);

        // Rotate
        #ifdef _USE_CUDA
        if(GPU)
            status = et_rotate_gpu(sourceImage, resultImage, angles[0],angles[1],angles[2], centers[0], centers[1], centers[2], background, axis_order);
        else
            status = et_rotate(sourceImage, resultImage, angles[0],angles[1],angles[2], centers[0], centers[1], centers[2], background, axis_order);
      	#else
      	    status = et_rotate(sourceImage, resultImage, angles[0],angles[1],angles[2], centers[0], centers[1], centers[2], background, axis_order);
        #endif
      	
	//Free
	if( sourceImage->fname != NULL ) free(sourceImage->fname) ;
	if( sourceImage->iname != NULL ) free(sourceImage->iname) ;
	(void)nifti_free_extensions( sourceImage ) ;
	free(sourceImage) ;
	
	if( resultImage->fname != NULL ) free(resultImage->fname) ;
	if( resultImage->iname != NULL ) free(resultImage->iname) ;
	(void)nifti_free_extensions( resultImage ) ;
	free(resultImage) ;

	return status;
}




extern "C" int et_array_convolve(float *image, int *image_size, float *out, int *out_size, float *psf, int *psf_size, int GPU)
{
        int status = 1;

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	//fprintf_verbose("\nS: %d %d %d",image_size[0],image_size[1],image_size[2]);
	nifti_image *imageImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        imageImage->data = (float *)(image);
	
	// Allocate the result nifti image
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];

        nifti_image *outImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        outImage->data = (float *)(out);

	//Allocate Point Spread Function nifti image
	//fprintf_verbose("\nPSF: %d %d %d",psf_size[0],psf_size[1],psf_size[2]);
	dim[1] = psf_size[0];
	dim[2] = psf_size[1];
	dim[3] = psf_size[2];
	nifti_image *psfImage;
        psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        psfImage->data = (float *)(psf);

        //Do convolution
        #ifdef _USE_CUDA
        if (GPU)
            status = et_convolve_gpu(imageImage, outImage, psfImage);
        else
            status = et_convolve(imageImage, outImage, psfImage);
        #else
            if (GPU)
                fprintf_verbose( "et_array_convolve: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_convolve(imageImage, outImage, psfImage);
        #endif

	//Free
	if( imageImage->fname != NULL ) free(imageImage->fname) ;
	if( imageImage->iname != NULL ) free(imageImage->iname) ;
	(void)nifti_free_extensions( imageImage ) ;
	free(imageImage) ;
	
	if( outImage->fname != NULL ) free(outImage->fname) ;
	if( outImage->iname != NULL ) free(outImage->iname) ;
	(void)nifti_free_extensions( outImage ) ;
	free(outImage) ;

	if( psfImage->fname != NULL ) free(psfImage->fname) ;
	if( psfImage->iname != NULL ) free(psfImage->iname) ;
	(void)nifti_free_extensions( psfImage ) ;
	free(psfImage) ;

	return status;
}



extern "C" int et_array_list_gpus(int *gpu_count, int *gpus_info_array)
{
        int status = 1;
	#ifdef _USE_CUDA
        status = et_list_gpus(gpu_count, gpus_info_array);
	#else
        gpu_count[0] = 0;
        status = niftyrec_error_nogpubuilt;
	#endif
        return status;
}


extern "C" int et_array_set_gpu_pointer(int* id)
{
        int status = 1;
        int ID = id[0]; 
        fprintf(stdout,"Requested GPU ID: %d\n",ID);
	#ifdef _USE_CUDA
        status = et_set_gpu(ID);
	#else
        status = niftyrec_error_nogpubuilt;
	#endif
        return status;
}


extern "C" int et_array_set_gpu(int id)
{
        int status = 1;
        fprintf(stdout,"Requested GPU ID: %d\n",id);
	#ifdef _USE_CUDA
        status = et_set_gpu(id);
	#else
        status = niftyrec_error_nogpubuilt;
	#endif
        return status;
}


extern "C" int et_array_isinstalled()
{
    return 0;
}


extern "C" int et_array_reset_gpu()
{
    #ifdef _USE_CUDA
        return et_reset_gpu();
    #else
      	return niftyrec_error_nogpubuilt;
    #endif
}


extern "C" int et_array_histogram_weighted(float *inputdata, float *weights, float *histogram, int N, int N_classes, int N_bins, double min_value, double max_value)
{
    int status = 1;
    float epsilon = 1e-8;
    int inputdata_int; 
    float sum  = (-min_value+5*epsilon); 
    float mult = (N_bins-10*epsilon)/(max_value-min_value); 
    float inputdata_scaled; 

    // clear histogram
    memset((void*) histogram, 0, N_bins*sizeof(float));

    // discretise and fill histogram 
    for (int i=0; i<N; i++)
        {
        inputdata_scaled = (inputdata[i]+sum)*mult; 
        if (inputdata_scaled<=epsilon) inputdata_scaled=epsilon;
        if (inputdata_scaled>=(N_bins-5*epsilon)) inputdata_scaled=N_bins-5*epsilon; 
        inputdata_int = floor(inputdata_scaled); 
        for (int k=0; k<N_classes; k++)
            {      
            //fprintf(stderr,"i:%d       k:%d     val:%f(%f-%d)     hist:%d       weight:%d\n",i,k,inputdata[i],inputdata_scaled,inputdata_int,k*N_bins+inputdata_int,k*N+i);
            histogram[k*N_bins+inputdata_int] = histogram[k*N_bins+inputdata_int] + weights[k*N+i]; 
            }
        }

    status = 0;
    return status; 
}



extern "C" unsigned short int et_array_is_block_multiple(unsigned short int size)
{
    return et_is_block_multiple(size); 
}



extern "C" unsigned short int et_array_get_block_size(void)
{
    return et_get_block_size();
}






//***************************************************************************************
/* SPECT */

extern "C" int SPECT_project_parallelholes(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float *Background, float *Background_attenuation, int *Gpu, int *Truncate_negative_values)
{
        int status = 0;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        int no_activity = 0;
        float *cameras_array;
        float *centers_array; 

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

        float background = *Background; 
        float background_attenuation = *Background_attenuation;
        int   GPU = *Gpu;
        int   truncate_negative_values = *Truncate_negative_values; 

        //activity or not? 
        if (activity_size[0] == 0 && activity_size[1] == 0 && activity_size[2] == 0) 
            {
            no_activity = 1;
            activity_size[0] = attenuation_size[0]; activity_size[1] = attenuation_size[1]; activity_size[2] = attenuation_size[2];
            }

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        if (no_attenuation && no_activity) 
            { 
            fprintf(stderr,"SPECT_project_parallelholes: Error - define at least one between 'activity' and 'attenuation'. \n"); 
            return niftyrec_error_parameters; 
            } 

        /* Check consistency of input */ 
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis. 
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3 || n_cameras_axis == 6))
            {
            fprintf_verbose("SPECT_project_parallelholes: Incorrect size of cameras %d %d. 'Cameras' must be [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6].\n",cameras_size[0],cameras_size[1]);
            return niftyrec_error_parameters;
            }

        //Size of psf must be odd and consistent with activity size
        if (!no_psf)
            {
            if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=activity_size[0])
                {
                fprintf_verbose("SPECT_project_parallelholes: 3D psf must be of size [h,k,N] for activity of size [N,m,N]; h,k odd.\n");
                return niftyrec_error_parameters;
                }
            }
/*
fprintf(stderr,"N_cameras:            %d \n",n_cameras);
fprintf(stderr,"N_cameras_axis:       %d \n",n_cameras_axis);
fprintf(stderr,"no_psf:               %d \n",no_psf);
fprintf(stderr,"no_attenuation:       %d \n",no_attenuation);
fprintf(stderr,"no_activity:          %d \n",no_activity);
fprintf(stderr,"background:           %f \n",background);
fprintf(stderr,"background_att:       %f \n",background_attenuation);
fprintf(stderr,"truncate:             %d \n",truncate_negative_values);
fprintf(stderr,"gpu:                  %d \n",GPU);
fprintf(stderr,"activity_size:        %d x %d x %d  \n",activity_size[0],activity_size[1],activity_size[2]);
fprintf(stderr,"sinogram_size:        %d x %d x %d  \n",sinogram_size[0],sinogram_size[1],sinogram_size[2]);
fprintf(stderr,"cameras_size:         %d x %d       \n",cameras_size[0],cameras_size[1] );
*/
        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if (cameras_array==NULL)
            return niftyrec_error_alloccpu;
        if (n_cameras_axis == 3 || n_cameras_axis == 6)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[1*n_cameras+cam] = cameras[cam];
            }

//        for(int i=0;i<n_cameras;i++)
//            fprintf(stderr,"Camera %d: %3.3f %3.3f %3.3f\n",i,cameras_array[0*n_cameras+i],cameras_array[1*n_cameras+i],cameras_array[2*n_cameras+i]);
            
        // Allocate array for center of rotation 
        centers_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if(n_cameras_axis == 6)
            {
            if (centers_array==NULL)
                return niftyrec_error_alloccpu;
            memcpy((void*) centers_array, (void*) &cameras[n_cameras*3], n_cameras*3*sizeof(float));
            }
        else
            {
            for (int cam=0; cam<n_cameras; cam++)
                {
                centers_array[0*n_cameras+cam] = (activity_size[0]-1)/2.0; 
                centers_array[1*n_cameras+cam] = (activity_size[1]-1)/2.0; 
                centers_array[2*n_cameras+cam] = (activity_size[2]-1)/2.0; 
                }
            }
//        for(int i=0;i<n_cameras;i++)
//            {
//            fprintf(stderr,"Camera %d: %3.3f %3.3f %3.3f\n",i,cameras_array[0*n_cameras+i],cameras_array[1*n_cameras+i],cameras_array[2*n_cameras+i]);
//            fprintf(stderr,"Center %d: %3.3f %3.3f %3.3f\n",i,centers_array[0*n_cameras+i],centers_array[1*n_cameras+i],centers_array[2*n_cameras+i]);  
//            }

	// Allocate nifti images
        int dim[8];
	dim[0]    = 3;
	dim[1]    = activity_size[0];
	dim[2]    = activity_size[1];
	dim[3]    = activity_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;

	// Allocate activity nifti images
	nifti_image *activityImage = NULL;
        if(!no_activity)
            {
            activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            activityImage->data = (float *)(activity);
            }

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result nifti image
	//fprintf_verbose( "\nN CAMERAS: %d ",n_cameras);

        dim[1] = activity_size[0];
        dim[2] = activity_size[1];
        dim[3] = n_cameras;	   

        nifti_image *sinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sinogramImage->data = (float *)(sinogram);

	//Allocate Point Spread Function nifti image
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float *)(psf);
            }

        //Do projection
        #ifdef _USE_CUDA
        if (GPU)
            status = et_project_gpu(activityImage, sinogramImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
        else
            status = et_project(activityImage, sinogramImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
        #else
            if (GPU)
                fprintf_verbose( "SPECT_project_parallelholes: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_project(activityImage, sinogramImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
        #endif

	//Free
        if(!no_activity)
            {
            if( activityImage->fname != NULL ) free(activityImage->fname) ;
            if( activityImage->iname != NULL ) free(activityImage->iname) ;
            (void)nifti_free_extensions( activityImage ) ;
            free(activityImage) ;
            }

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }
	
	if( sinogramImage->fname != NULL ) free(sinogramImage->fname) ;
	if( sinogramImage->iname != NULL ) free(sinogramImage->iname) ;
	(void)nifti_free_extensions( sinogramImage ) ;
	free(sinogramImage) ;

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);
        free(centers_array);

	return status;
}





extern "C" int SPECT_backproject_parallelholes(float *sino, int *sino_size, float *bkpr, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float *Background, float *Background_attenuation, int *Gpu, int *Truncate_negative_values)
{
	int status;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;
        float *centers_array; 
        int activity_size[3]; 

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

        float background = *Background; 
        float background_attenuation = *Background_attenuation;
        int   GPU = *Gpu;
        int   truncate_negative_values = *Truncate_negative_values; 

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3 || n_cameras_axis == 6))
            {
            fprintf_verbose("SPECT_backproject_parallelholes: 'Cameras' must be [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6]\n");
            return status;
            }

        //Sino must be of size [Nxmxn_cameras]; 
        if (sino_size[2] != n_cameras)
            {
            fprintf_verbose("SPECT_backproject_parallelholes: 3D sino must be of size [N,m,n_cameras].\n");
            return status;
            }
        //Size of psf must be odd and consistent with activity size
        if (!no_psf)
            {
            if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=sino_size[0])
                {
                fprintf_verbose("SPECT_backproject_parallelholes: 3D psf must be of size [h,k,N] for activity of size [N,m,N]; h,k odd.\n");
                return status;
                }
            }

        activity_size[0] = sino_size[0]; 
        activity_size[1] = sino_size[1]; 
        activity_size[2] = sino_size[0]; 

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (cameras_array==NULL)
            return niftyrec_error_alloccpu; 
        if (n_cameras_axis == 3 || n_cameras_axis == 6)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[1*n_cameras+cam] = cameras[cam];
            }

        // Allocate array for center of rotation 
        centers_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if(n_cameras_axis == 6)
            {
            if (centers_array==NULL)
                return niftyrec_error_alloccpu;
            memcpy((void*) centers_array, (void*) &cameras[n_cameras*3], n_cameras*3*sizeof(float));
            }
        else
            {
            for (int cam=0; cam<n_cameras; cam++)
                {
                centers_array[0*n_cameras+cam] = (activity_size[0]-1)/2.0; 
                centers_array[1*n_cameras+cam] = (activity_size[1]-1)/2.0; 
                centers_array[2*n_cameras+cam] = (activity_size[2]-1)/2.0; 
                }
            }

	// Allocate backprojection (result) image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = bkpr_size[0];
	dim[2]    = bkpr_size[1];
	dim[3]    = bkpr_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *bkprImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        bkprImage->data = (float*) bkpr;
	
        // Allocate attenuation image 
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float*) attenuation;
            }

	// Allocate the sinogram (input) image
	dim[3]    = n_cameras;
        nifti_image *sinoImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sinoImage->data = (float*) sino;

	//Allocate Point Spread Function
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float*) psf;
            }
	//Backproject
	#ifdef _USE_CUDA
	if(GPU)
	    status = et_backproject_gpu(sinoImage, bkprImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
	else
	    status = et_backproject(sinoImage, bkprImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);	    
	#else
	    status = et_backproject(sinoImage, bkprImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);	    
	#endif
	
	//Free (free nifti images but not their data arrays)
	if( sinoImage->fname != NULL ) free(sinoImage->fname) ;
	if( sinoImage->iname != NULL ) free(sinoImage->iname) ;
	(void)nifti_free_extensions( sinoImage ) ;
	free(sinoImage) ;
	
	if( bkprImage->fname != NULL ) free(bkprImage->fname) ;
	if( bkprImage->iname != NULL ) free(bkprImage->iname) ;
	(void)nifti_free_extensions( bkprImage ) ;
	free(bkprImage) ;

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);
        free(centers_array); 

	return status;
}



extern "C" int et_array_gradient_attenuation(float *sino, int *sino_size, float *activity, int *activity_size, float *gradient, int *gradient_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU, int truncate_negative_values)
{
	int status;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        float *cameras_array;
        float *centers_array; 

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (sino_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3 || n_cameras_axis == 6))
            {
            fprintf_verbose("et_array_gradient_attenuation: 'Cameras' must be [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6]\n");
            return status;
            }
        if (dims==2)
            //Sino must be of size [Nxn_cameras]
            {
            if (sino_size[1] != n_cameras)
                {
                fprintf_verbose("et_array_gradient_attenuation: 2D sinogram must be of size [N,n_cameras].\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]!=sino_size[0])
                    {
                    fprintf_verbose("et_array_gradient_attenuation: 2D psf must be of size [h,N]; h odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Sino must be of size [Nxmxn_cameras]; 
            {
            if (sino_size[2] != n_cameras)
                {
                fprintf_verbose("et_array_gradient_attenuation: 3D sino must be of size [N,m,n_cameras].\n");
                return status;
                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=sino_size[0])
                    {
                    fprintf_verbose("et_array_gradient_attenuation: 3D psf must be of size [h,k,N] for activity of size [N,m,N]; h,k odd.\n");
                    return status;
                    }
                }
            }

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (cameras_array==NULL)
            return niftyrec_error_alloccpu; 
        if (n_cameras_axis == 3 || n_cameras_axis == 6)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[1*n_cameras+cam] = cameras[cam];
            }

        // Allocate array for center of rotation 
        centers_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if(n_cameras_axis == 6)
            {
            if (centers_array==NULL)
                return niftyrec_error_alloccpu;
            memcpy((void*) centers_array, (void*) &cameras[n_cameras*3], n_cameras*3*sizeof(float));
            }
        else
            {
            for (int cam=0; cam<n_cameras; cam++)
                {
                centers_array[0*n_cameras+cam] = (activity_size[0]-1)/2.0; 
                centers_array[1*n_cameras+cam] = (activity_size[1]-1)/2.0; 
                centers_array[2*n_cameras+cam] = (activity_size[2]-1)/2.0; 
                }
            }

	// Allocate backprojection (result) image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = gradient_size[0];
	dim[2]    = gradient_size[1];
	dim[3]    = gradient_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *gradientImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        gradientImage->data = (float*) gradient;
	
        // Allocate attenuation image 
        nifti_image *attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        attenuationImage->data = (float*) attenuation;

        // Allocate activity image 
	nifti_image *activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        activityImage->data = (float *)(activity);

	// Allocate the sinogram (input) image
	dim[1]    = gradient_size[0];
	dim[2]    = gradient_size[1];
	dim[3]    = n_cameras;
        nifti_image *sinoImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sinoImage->data = (float*) sino;

	//Allocate Point Spread Function
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float*) psf;
            }
	//Backproject
	#ifdef _USE_CUDA
	if(GPU)
	    status = et_gradient_attenuation_gpu(gradientImage, sinoImage, activityImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
	else
	    status = et_gradient_attenuation(gradientImage, sinoImage, activityImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
	#else
	    status = et_gradient_attenuation(gradientImage, sinoImage, activityImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation, truncate_negative_values);
	#endif
	
	//Free (free nifti images but not their data arrays)
	if( sinoImage->fname != NULL ) free(sinoImage->fname) ;
	if( sinoImage->iname != NULL ) free(sinoImage->iname) ;
	(void)nifti_free_extensions( sinoImage ) ;
	free(sinoImage) ;
	
	if( gradientImage->fname != NULL ) free(gradientImage->fname) ;
	if( gradientImage->iname != NULL ) free(gradientImage->iname) ;
	(void)nifti_free_extensions( gradientImage ) ;
	free(gradientImage) ;

        if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
        if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
        (void)nifti_free_extensions( attenuationImage ) ;
        free(attenuationImage) ;

        if( activityImage->fname != NULL ) free(activityImage->fname) ;
        if( activityImage->iname != NULL ) free(activityImage->iname) ;
        (void)nifti_free_extensions( activityImage ) ;
        free(activityImage) ;

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);
        free(centers_array);

	return status;
}



extern "C" int et_array_calculate_size_psf(unsigned int *psf_size_x, unsigned int *psf_size_y, float fwhm_pixels_dist0, float sensitivity0, float dist0, float fwhm_pixels_dist1, float sensitivity1, float dist1)
{
    *psf_size_x = 5;
    *psf_size_y = 5;
    return 0;
}


extern "C" int et_array_make_psf(float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float fwhm_pixels_dist0, float sensitivity0, float dist0, float fwhm_pixels_dist1, float sensitivity1, float dist1, unsigned int N_psf_planes)
{
    int status;
    unsigned int psf_size_x_calc;
    unsigned int psf_size_y_calc;
    status = et_array_calculate_size_psf(&psf_size_x_calc, &psf_size_y_calc, fwhm_pixels_dist0, sensitivity0, dist0, fwhm_pixels_dist1, sensitivity1, dist1); 
    if (status!=0 || psf_size_x != psf_size_x_calc || psf_size_y != psf_size_y_calc)
        {
        fprintf(stderr,"et_array_make_psf: PSF size mismatch. \n");
        return 1;
        }

    for (int i=0; i<psf_size_x*psf_size_y*N_psf_planes; i++)
        psf_data[i] = 1;

    return 0; 
}


extern "C" int et_array_make_cameras(float *cameras_data, float firstcamera_deg, float lastcamera_deg, unsigned int n_cameras, unsigned int rotation_axis)
{
    for (int cam=0; cam<n_cameras; cam++)
        {
        float angle = (firstcamera_deg + (lastcamera_deg-firstcamera_deg)*cam / (n_cameras-1))/180.0f*PI; 
        for (unsigned int axis=0; axis<3; axis++)
            {
            if (axis==rotation_axis)
                cameras_data[cam+axis*n_cameras] = angle;
            else
                cameras_data[cam+axis*n_cameras] = 0.0f; 
            }
        }
    return 0;
}


extern "C" int et_array_osem(float *activity_data, unsigned int size_x, unsigned int size_y, unsigned int subset_order, float *sinogram_data, int n_cameras, float firstcamera, float lastcamera, unsigned int rotation_axis, int iterations, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu)
{
    int status = 0;

    for (int i=0; i<size_x*size_y*size_x; i++)
        activity_data[i] = 1.0f; 

    float *cameras_data = (float*) malloc(n_cameras*3*sizeof(float)); 
    if (cameras_data==NULL)
        return niftyrec_error_alloccpu; 
    status = status + et_array_make_cameras(cameras_data, firstcamera, lastcamera, n_cameras, rotation_axis); 

    for (int iter_osem=0; iter_osem<iterations; iter_osem++)
        {
        fprintf(stderr,"et_array_osem: step %d/%d \n",iter_osem,iterations); 
        status = status + et_array_osem_step(activity_data, size_x, size_y, subset_order, sinogram_data, n_cameras, cameras_data, use_attenuation, attenuation_data, use_psf, psf_data, psf_size_x, psf_size_y, background, background_attenuation, epsilon, use_gpu); 
        }

    free(cameras_data); 
    return status; 
}

float random_0_1()
{
  float scale=RAND_MAX+1.;
  float base=rand()/scale;
  float fine=rand()/scale;
  return base+fine/scale;
}



extern "C" int et_array_osem_step(float *activity_data, int size_x, int size_y, unsigned int subset_order, float *sinogram_data, int n_cameras, float *cameras_data, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu) 
{
    int status=0;
    int partial_status; 
    int sinogram_size[3]; sinogram_size[0]=size_x; sinogram_size[1]=size_y; sinogram_size[2]=n_cameras;
    int activity_size[3]; activity_size[0]=sinogram_size[0]; activity_size[1]=sinogram_size[1]; activity_size[2]=sinogram_size[0]; 
    int psf_size[3] = {0,0,0}; if (use_psf) {psf_size[0]=psf_size_x; psf_size[1]=psf_size_y; psf_size[2]=size_x; }; 
    int attenuation_size[3] = {0,0,0}; if (use_attenuation) {attenuation_size[0]=activity_size[0]; attenuation_size[1]=activity_size[1]; attenuation_size[2]=activity_size[2]; }; 

    // Select random subset: select randomly the first camera, 
    // then select successive camera by drawing from a Gaussian 
    // centred 'subset_order' indexes away 
    // and procede that way in the same direction until 'N_cameras/subset_order' 
    // cameras are selected

    unsigned int n_cameras_sub = n_cameras / subset_order; 
    unsigned int *cameras_indexes = (unsigned int*) malloc(n_cameras_sub*sizeof(unsigned int)); 
    int sinogram_sub_size[3]; sinogram_sub_size[0]=size_x; sinogram_sub_size[1]=size_y; sinogram_sub_size[2]=n_cameras_sub;
    int cameras_sub_size[2]; cameras_sub_size[1]=3; cameras_sub_size[0]=n_cameras_sub; 
    float *projection_sub_data = (float*) malloc(sinogram_sub_size[0]*sinogram_sub_size[1]*sinogram_sub_size[2]*sizeof(float)); 
    float *poisson_gradient_data = (float*) malloc(activity_size[0]*activity_size[1]*activity_size[2]*sizeof(float));
    float *normalisation_data = (float*) malloc(activity_size[0]*activity_size[1]*activity_size[2]*sizeof(float));  
    float *sinogram_sub_ones_data = (float*) malloc(sinogram_sub_size[0]*sinogram_sub_size[1]*sinogram_sub_size[2]*sizeof(float)); 
    float *cameras_sub_data = (float*) malloc(n_cameras_sub*3*sizeof(float));
    float *sinogram_sub_data = (float*) malloc(sinogram_sub_size[0]*sinogram_sub_size[1]*sinogram_sub_size[2]*sizeof(float)); 

    for (int i=0; i<sinogram_sub_size[0]*sinogram_sub_size[1]*sinogram_sub_size[2]; i++) 
        sinogram_sub_ones_data[i]=1.0f; 

    for (int cam=0; cam<n_cameras_sub ; cam++)
        {
        float R = (float)random_0_1();
        cameras_indexes[cam] = n_cameras*R; 
        fprintf(stderr,"R: %f  Camera index: %d\n",R,cameras_indexes[cam]);
        if(cameras_indexes[cam]>=n_cameras) cameras_indexes[cam]=n_cameras;
        for (unsigned int axis=0; axis<3; axis++)
             cameras_sub_data[cam+axis*n_cameras_sub] = cameras_data[cameras_indexes[cam]+axis*n_cameras]; 
        }

    for (int cam=0; cam<n_cameras_sub; cam++)
        {
        for (int i=0; i<sinogram_size[0]*sinogram_size[1]; i++)
            sinogram_sub_data[cam*sinogram_size[0]*sinogram_size[1]+i]=sinogram_data[cameras_indexes[cam]*sinogram_size[0]*sinogram_size[1]+i];
        }

    int truncate_negative_values=1;
    partial_status = SPECT_project_parallelholes(activity_data, activity_size, projection_sub_data, sinogram_sub_size, cameras_sub_data, cameras_sub_size, psf_data, psf_size, attenuation_data, attenuation_size, &background, &background_attenuation, &use_gpu, &truncate_negative_values);     
    status = status + partial_status; 
    if (partial_status)
        fprintf(stderr, "et_array_osem_step: error while performing projection \n");
    for (int i=0; i<sinogram_sub_size[0]*sinogram_sub_size[1]*sinogram_sub_size[2]; i++) 
        {
        if(projection_sub_data[i]<epsilon)
            projection_sub_data[i]=sinogram_sub_data[i]/epsilon; 
        else
            projection_sub_data[i]=sinogram_sub_data[i]/projection_sub_data[i];
        }

    partial_status = SPECT_backproject_parallelholes(projection_sub_data, sinogram_sub_size, poisson_gradient_data, activity_size, cameras_sub_data, cameras_sub_size, psf_data, psf_size, attenuation_data, attenuation_size, &background, &background_attenuation, &use_gpu, &truncate_negative_values); 
    status = status + partial_status; 
    if (partial_status)
        fprintf(stderr, "et_array_osem_step: error while performing back-projection \n");

    partial_status = SPECT_backproject_parallelholes(sinogram_sub_ones_data, sinogram_sub_size, normalisation_data, activity_size, cameras_sub_data, cameras_sub_size, psf_data, psf_size, attenuation_data, attenuation_size, &background, &background_attenuation, &use_gpu, &truncate_negative_values); 
    status = status + partial_status; 
    if (partial_status)
        fprintf(stderr, "et_array_osem_step: error while computing normalisation \n");
    for (int i=0; i<activity_size[0]*activity_size[1]*activity_size[2]; i++)
        if(normalisation_data[i]<epsilon)
            normalisation_data[i]=epsilon; 

    for (int i=0; i<activity_size[0]*activity_size[1]*activity_size[2]; i++)
        activity_data[i] *= poisson_gradient_data[i] / normalisation_data[i]; 

    free(projection_sub_data);
    free(poisson_gradient_data); 
    free(normalisation_data); 
    free(sinogram_sub_ones_data); 
    free(cameras_indexes);
    free(cameras_sub_data); 
    free(sinogram_sub_data);
    return status; 
}



extern "C" int et_array_mlem(float *activity_data, unsigned int size_x, unsigned int size_y, float *sinogram_data, int n_cameras, float firstcamera, float lastcamera, unsigned int rotation_axis, int iterations, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu)
{
    int status = 0;

    for (int i=0; i<size_x*size_y*size_x; i++)
        activity_data[i] = 1.0f; 

    float *cameras_data = (float*) malloc(n_cameras*3*sizeof(float)); 
    status = status + et_array_make_cameras(cameras_data, firstcamera, lastcamera, n_cameras, rotation_axis); 

    for (int iter_mlem=0; iter_mlem<iterations; iter_mlem++)
        {
        fprintf(stderr,"et_array_mlem: step %d/%d \n",iter_mlem,iterations); 
        status = status + et_array_mlem_step(activity_data, size_x, size_y, sinogram_data, n_cameras, cameras_data, use_attenuation, attenuation_data, use_psf, psf_data, psf_size_x, psf_size_y, background, background_attenuation, epsilon, use_gpu); 
        }

    free(cameras_data); 
    return status; 
}



extern "C" int et_array_mlem_step(float *activity_data, int size_x, int size_y, float *sinogram_data, int n_cameras, float *cameras_data, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu) 
{
    int status=0;
    int partial_status; 
    int sinogram_size[3]; sinogram_size[0]=size_x; sinogram_size[1]=size_y; sinogram_size[2]=n_cameras;
    int activity_size[3]; activity_size[0]=sinogram_size[0]; activity_size[1]=sinogram_size[1]; activity_size[2]=sinogram_size[0]; 
    int cameras_size[2]; cameras_size[1]=3; cameras_size[0]=n_cameras; 
    int psf_size[3] = {0,0,0}; if (use_psf) {psf_size[0]=psf_size_x; psf_size[1]=psf_size_y; psf_size[2]=size_x; }; 
    int attenuation_size[3] = {0,0,0}; if (use_attenuation) {attenuation_size[0]=activity_size[0]; attenuation_size[1]=activity_size[1]; attenuation_size[2]=activity_size[2]; }; 

    float *projection_data = (float*) malloc(sinogram_size[0]*sinogram_size[1]*sinogram_size[2]*sizeof(float)); 
    float *poisson_gradient_data = (float*) malloc(activity_size[0]*activity_size[1]*activity_size[2]*sizeof(float));
    float *normalisation_data = (float*) malloc(activity_size[0]*activity_size[1]*activity_size[2]*sizeof(float));  
    float *sinogram_ones_data = (float*) malloc(sinogram_size[0]*sinogram_size[1]*sinogram_size[2]*sizeof(float)); 
    for (int i=0; i<sinogram_size[0]*sinogram_size[1]*sinogram_size[2]; i++) 
        sinogram_ones_data[i]=1.0f; 
    
    int truncate_negative_values=1;
    partial_status = SPECT_project_parallelholes(activity_data, activity_size, projection_data, sinogram_size, cameras_data, cameras_size, psf_data, psf_size, attenuation_data, attenuation_size, &background, &background_attenuation, &use_gpu, &truncate_negative_values);     
    status = status + partial_status; 
    if (partial_status)
        fprintf(stderr, "et_array_mlem_step: error while performing projection \n");
    for (int i=0; i<sinogram_size[0]*sinogram_size[1]*sinogram_size[2]; i++) 
        {
        if(projection_data[i]<epsilon)
            projection_data[i]=sinogram_data[i]/epsilon; 
        else
            projection_data[i]=sinogram_data[i]/projection_data[i];
        }

    partial_status = SPECT_backproject_parallelholes(projection_data, sinogram_size, poisson_gradient_data, activity_size, cameras_data, cameras_size, psf_data, psf_size, attenuation_data, attenuation_size, &background, &background_attenuation, &use_gpu, &truncate_negative_values); 
    status = status + partial_status; 
    if (partial_status)
        fprintf(stderr, "et_array_mlem_step: error while performing back-projection \n");

    partial_status = SPECT_backproject_parallelholes(sinogram_ones_data, sinogram_size, normalisation_data, activity_size, cameras_data, cameras_size, psf_data, psf_size, attenuation_data, attenuation_size, &background, &background_attenuation, &use_gpu, &truncate_negative_values); 
    status = status + partial_status; 
    if (partial_status)
        fprintf(stderr, "et_array_mlem_step: error while computing normalisation \n");
    for (int i=0; i<activity_size[0]*activity_size[1]*activity_size[2]; i++)
        if(normalisation_data[i]<epsilon)
            normalisation_data[i]=epsilon; 

    for (int i=0; i<activity_size[0]*activity_size[1]*activity_size[2]; i++)
        activity_data[i] *= poisson_gradient_data[i] / normalisation_data[i]; 

    free(projection_data);
    free(poisson_gradient_data); 
    free(normalisation_data); 
    free(sinogram_ones_data); 
    return status; 
}



extern "C" int et_array_fisher_grid(float *activity_ptr, int *activity_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU)
{
        int from_projection = 0;
	int status = 1;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;
        float *centers_array;

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (activity_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3 || n_cameras_axis == 6))
            {
            fprintf(stderr,"et_array_fisher_grid: Incorrect size of cameras %d %d. 'Cameras' must be [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6].\n",cameras_size[0],cameras_size[1]);
            return status;
            }
        if (dims==2)
            //Activity must be of size [NxN]
            {
            if (activity_size[0] != activity_size[1])
                {
                fprintf(stderr,"et_array_fisher_grid: 2D activity must be of size [N,N].\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1)
                    {
                    fprintf(stderr,"et_array_fisher_grid: 2D psf must be of size [h,k]; h,k odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Activity must be of size [NxmxN]; 
            {
            if (activity_size[0] != activity_size[2])
                {
                fprintf(stderr,"et_array_fisher_grid: 3D activity must be of size [N,m,N]\n");
                return status;
                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=activity_size[0])
                    {
                    fprintf(stderr,"et_array_fisher_grid: 3D psf must be of size [h,k,N] for activity of size [N,m,N]; h,k odd.\n");
                    return status;
                    }
                }
            }

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (cameras_array==NULL)
            return niftyrec_error_alloccpu; 
        if (n_cameras_axis == 3 || n_cameras_axis == 6)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[1*n_cameras+cam] = cameras[cam];
            }

        // Allocate array for center of rotation 
        centers_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if(n_cameras_axis == 6)
            {
            if (centers_array==NULL)
                return niftyrec_error_alloccpu;
            memcpy((void*) centers_array, (void*) &cameras[n_cameras*3], n_cameras*3*sizeof(float));
            }
        else
            {
            for (int cam=0; cam<n_cameras; cam++)
                {
                centers_array[0*n_cameras+cam] = (activity_size[0]-1)/2.0; 
                centers_array[1*n_cameras+cam] = (activity_size[1]-1)/2.0; 
                centers_array[2*n_cameras+cam] = (activity_size[2]-1)/2.0; 
                }
            }

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = activity_size[0];
	dim[2]    = activity_size[1];
	dim[3]    = activity_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	//fprintf(stderr,"\nS: %d %d %d",activity_size[0],activity_size[1],activity_size[2]);
	nifti_image *activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        activityImage->data = (float *)(activity_ptr);

        // Allocate grid nifti image
        nifti_image *gridImage = NULL;
        gridImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        gridImage->data = (float *)(grid_ptr);        

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result Fisher nifti image
        dim[1] = fisher_size[0];
        dim[2] = fisher_size[1];
        dim[3] = 1;
        nifti_image *fisherImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        fisherImage->data = (float *)(fisher_ptr);

        nifti_image *fisherpriorImage = NULL;
        if (fisher_prior_ptr!=NULL)
            {
            fisherpriorImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            fisherpriorImage->data = (float *)(fisher_prior_ptr);
            }

	//Allocate Point Spread Function nifti image
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float *)(psf);
            }

        //Compute Fisher Information Matrix
        #ifdef _USE_CUDA
        if (GPU)
            status = et_fisher_grid_gpu(from_projection, activityImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation); 
        else
            status = et_fisher_grid(from_projection, activityImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation); 
        #else
            if (GPU)
                fprintf(stderr, "et_array_fisher_grid: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_fisher_grid(from_projection, activityImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation); 
        #endif

	//Free
	if( activityImage->fname != NULL ) free(activityImage->fname) ;
	if( activityImage->iname != NULL ) free(activityImage->iname) ;
	(void)nifti_free_extensions( activityImage ) ;
	free(activityImage) ;

	(void)nifti_free_extensions( gridImage ) ;
	free(gridImage) ;

	(void)nifti_free_extensions( fisherImage ) ;
	free(fisherImage) ;

        if (fisher_prior_ptr!=NULL)
            {         
            (void)nifti_free_extensions( fisherpriorImage ) ;
            free(fisherpriorImage) ;
            }

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);
        free(centers_array);

	return status;
}



extern "C" int et_array_fisher_grid_projection(float *sinogram_ptr, int *sinogram_size, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU)
{
        int from_projection = 1;
	int status = 1;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;
        float *centers_array; 
        int activity_size[3]; 

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (bkpr_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3 || n_cameras_axis == 6))
            {
            fprintf(stderr,"et_array_fisher_grid_projection: Incorrect size of cameras %d %d. 'Cameras' must be [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6].\n",cameras_size[0],cameras_size[1]);
            return status;
            }
        if (dims==2)
            //Activity must be of size [NxN]
            {
            if (bkpr_size[0] != bkpr_size[1])
                {
                fprintf(stderr,"et_array_fisher_grid_projection: 2D activity must be of size [N,N].\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1)
                    {
                    fprintf(stderr,"et_array_fisher_grid_projection: 2D psf must be of size [h,k]; h,k odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Activity must be of size [NxmxN]
            {
//            if (bkpr_size[0] != bkpr_size[1])
//                {
//                fprintf(stderr,"et_array_fisher_grid_projection: 3D activity must be of size [N,m,N].\n");
//                return status;
//                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=bkpr_size[0])
                    {
                    fprintf(stderr,"et_array_fisher_grid_projection: 3D psf must be of size [h,k,N] for activity of size [N,m,N]; h,k odd.\n");
                    return status;
                    }
                }
            }

        activity_size[0]=sinogram_size[0];
        activity_size[1]=sinogram_size[1];
        activity_size[2]=sinogram_size[0];

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (cameras_array==NULL)
            return niftyrec_error_alloccpu; 
        if (n_cameras_axis == 3 || n_cameras_axis == 6)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[1*n_cameras+cam] = cameras[cam];
            }

        // Allocate array for center of rotation 
        centers_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if(n_cameras_axis == 6)
            {
            if (centers_array==NULL)
                return niftyrec_error_alloccpu;
            memcpy((void*) centers_array, (void*) &cameras[n_cameras*3], n_cameras*3*sizeof(float));
            }
        else
            {
            for (int cam=0; cam<n_cameras; cam++)
                {
                centers_array[0*n_cameras+cam] = (activity_size[0]-1)/2.0; 
                centers_array[1*n_cameras+cam] = (activity_size[1]-1)/2.0; 
                centers_array[2*n_cameras+cam] = (activity_size[2]-1)/2.0; 
                }
            }

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = sinogram_size[0];
	dim[2]    = sinogram_size[1];
	dim[3]    = sinogram_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *projectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        projectionImage->data = (float *)(sinogram_ptr);

        // Allocate grid nifti image
	dim[1]    = bkpr_size[0];
	dim[2]    = bkpr_size[1];
	dim[3]    = bkpr_size[2];
        nifti_image *gridImage = NULL;
        gridImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        gridImage->data = (float *)(grid_ptr);        

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result Fisher nifti image
        dim[1] = fisher_size[0];
        dim[2] = fisher_size[1];
        dim[3] = 1;
        nifti_image *fisherImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        fisherImage->data = (float *)(fisher_ptr);

        nifti_image *fisherpriorImage = NULL;
        if (fisher_prior_ptr!=NULL)
            {
            fisherpriorImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            fisherpriorImage->data = (float *)(fisher_prior_ptr);
            }

	//Allocate Point Spread Function nifti image
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float *)(psf);
            }

        //Compute Fisher Information Matrix
        #ifdef _USE_CUDA
        if (GPU)
            status = et_fisher_grid_gpu(from_projection, projectionImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation); 
        else
            status = et_fisher_grid(from_projection, projectionImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation); 
        #else
            if (GPU)
                fprintf(stderr,"et_array_fisher_grid_projection: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_fisher_grid(from_projection, projectionImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, background_attenuation); 
        #endif

	//Free
	if( projectionImage->fname != NULL ) free(projectionImage->fname) ;
	if( projectionImage->iname != NULL ) free(projectionImage->iname) ;
	(void)nifti_free_extensions( projectionImage ) ;
	free(projectionImage) ;

	(void)nifti_free_extensions( gridImage ) ;
	free(gridImage) ;

	(void)nifti_free_extensions( fisherImage ) ;
	free(fisherImage) ;

        if (fisher_prior_ptr!=NULL)
            {         
            (void)nifti_free_extensions( fisherpriorImage ) ;
            free(fisherpriorImage) ;
            }

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);
        free(centers_array);

	return status;
}




extern "C" int et_array_project_partial(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *partialsum, int *partialsum_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float *background_image, float background_attenuation, int GPU, int truncate_negative_values, int do_rotate_partial)
{
        int status = 0;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        int no_activity = 0;
        int no_background_image = 0;
        float *cameras_array;
        float *centers_array; 

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

        //activity or not? 
        if (activity_size[0] == 0 && activity_size[1] == 0 && activity_size[2] == 0) 
            {
            no_activity = 1;
            activity_size[0] = attenuation_size[0]; activity_size[1] = attenuation_size[1]; activity_size[2] = attenuation_size[2];
            }

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        if (no_attenuation && no_activity) 
            { 
            fprintf(stderr,"et_array_project_partial: Error - define at least one between 'activity' and 'attenuation'. \n"); 
            return niftyrec_error_parameters; 
            } 

        //background image or not?
        no_background_image = 0; 
        if (background_image==NULL)
            no_background_image = 1; 

        /* Check consistency of input */ 
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis. 
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3 || n_cameras_axis == 6))
            {
            fprintf_verbose("et_array_project_partial: Incorrect size of cameras %d %d. 'Cameras' must be [n_cameras x 1] or [n_cameras x 3] or [n_cameras x 6].\n",cameras_size[0],cameras_size[1]);
            return niftyrec_error_parameters;
            }

        //Size of psf must be odd and consistent with activity size
        if (!no_psf)
            {
            if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=activity_size[0])
                {
                fprintf_verbose("et_array_project_partial: 3D psf must be of size [h,k,N] for activity of size [N,m,N]; h,k odd.\n");
                return niftyrec_error_parameters;
                }
            }


        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if (cameras_array==NULL)
            return niftyrec_error_alloccpu;
        if (n_cameras_axis == 3 || n_cameras_axis == 6)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[1*n_cameras+cam] = cameras[cam];
            }

        // Allocate array for center of rotation 
        centers_array = (float *)malloc(n_cameras*3*sizeof(float)); 
        if(n_cameras_axis == 6)
            {
            if (centers_array==NULL)
                return niftyrec_error_alloccpu;
            memcpy((void*) centers_array, (void*) &cameras[n_cameras*3], n_cameras*3*sizeof(float));
            }
        else
            {
            for (int cam=0; cam<n_cameras; cam++)
                {
                centers_array[0*n_cameras+cam] = (activity_size[0]-1)/2.0; 
                centers_array[1*n_cameras+cam] = (activity_size[1]-1)/2.0; 
                centers_array[2*n_cameras+cam] = (activity_size[2]-1)/2.0; 
                }
            }

	// Allocate nifti images
        int dim[8];
	dim[0]    = 3;
	dim[1]    = activity_size[0];
	dim[2]    = activity_size[1];
	dim[3]    = activity_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;

	// Allocate activity nifti images
	nifti_image *activityImage = NULL;
        if(!no_activity)
            {
            activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            activityImage->data = (float *)(activity);
            }

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result nifti image
	//fprintf_verbose( "\nN CAMERAS: %d ",n_cameras);
        dim[1] = activity_size[0];
        dim[2] = activity_size[1];
        dim[3] = n_cameras;	   

        nifti_image *sinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sinogramImage->data = (float *)(sinogram);

	//Allocate Point Spread Function nifti image
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float *)(psf);
            }

	// Allocate the nifti image for the partial sum
        dim[0] = 4;
        dim[1] = partialsum_size[0];
        dim[2] = partialsum_size[1];
        dim[3] = partialsum_size[2];
        dim[4] = partialsum_size[3];

        nifti_image *partialsumImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        partialsumImage->data = (float *)(partialsum);

        // Allocate background nifti image
        nifti_image *backgroundImage = NULL;
        if(!no_background_image)
            {
            dim[0]    = 2;
            dim[1]    = activity_size[0];
            dim[2]    = activity_size[1];
            dim[3]    = 1;
            dim[4]    = 1;
            backgroundImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            backgroundImage->data = (float *)(background_image);        
            }

        //Do projection
        #ifdef _USE_CUDA
        if (GPU)
            status = et_project_partial_gpu(activityImage, sinogramImage, partialsumImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, backgroundImage, background_attenuation, truncate_negative_values, do_rotate_partial);
        else
            status = et_project_partial(activityImage, sinogramImage, partialsumImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, backgroundImage, background_attenuation, truncate_negative_values, do_rotate_partial);
        #else
            if (GPU)
                fprintf_verbose( "et_array_project_partial: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_project_partial(activityImage, sinogramImage, partialsumImage, psfImage, attenuationImage, cameras_array, centers_array, n_cameras, background, backgroundImage, background_attenuation, truncate_negative_values, do_rotate_partial);
        #endif

	//Free
        if(!no_activity)
            {
            if( activityImage->fname != NULL ) free(activityImage->fname) ;
            if( activityImage->iname != NULL ) free(activityImage->iname) ;
            (void)nifti_free_extensions( activityImage ) ;
            free(activityImage) ;
            }

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }
	
        if(!no_background_image)
            {
            if( backgroundImage->fname != NULL ) free(backgroundImage->fname) ;
            if( backgroundImage->iname != NULL ) free(backgroundImage->iname) ;
            (void)nifti_free_extensions( backgroundImage ) ;
            free(backgroundImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

	if( sinogramImage->fname != NULL ) free(sinogramImage->fname) ;
	if( sinogramImage->iname != NULL ) free(sinogramImage->iname) ;
	(void)nifti_free_extensions( sinogramImage ) ;
	free(sinogramImage) ;

	if( partialsumImage->fname != NULL ) free(partialsumImage->fname) ;
	if( partialsumImage->iname != NULL ) free(partialsumImage->iname) ;
	(void)nifti_free_extensions( partialsumImage ) ;
	free(partialsumImage) ;

        free(cameras_array);
        free(centers_array);

	return status;
}




//***************************************************************************************
/* PET */

extern "C" int PET_project(int *param)
{
    return STATUS_SUCCESS; 
}

extern "C" unsigned int PET_backproject(void)
{
    return STATUS_SUCCESS; 
}

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
                           unsigned int *direction, unsigned int *block_size, float *time_profiling)
{
    fprintf(stderr,"Compressed projection ..\n");
if (0) {    
    fprintf(stderr,"  - N_activity_x:                 %d \n",*N_activity_x);
    fprintf(stderr,"  - N_activity_y:                 %d \n",*N_activity_y);
    fprintf(stderr,"  - N_activity_z:                 %d \n",*N_activity_z);
    fprintf(stderr,"  - activity_size_x:              %f \n",*activity_size_x);
    fprintf(stderr,"  - activity_size_y:              %f \n",*activity_size_y);
    fprintf(stderr,"  - activity_size_z:              %f \n",*activity_size_z);
    fprintf(stderr,"  - T_activity_x:                 %f \n",*T_activity_x); 
    fprintf(stderr,"  - T_activity_y:                 %f \n",*T_activity_y); 
    fprintf(stderr,"  - T_activity_z:                 %f \n",*T_activity_z); 
    fprintf(stderr,"  - R_activity_x:                 %f \n",*R_activity_x); 
    fprintf(stderr,"  - R_activity_y:                 %f \n",*R_activity_y); 
    fprintf(stderr,"  - R_activity_z:                 %f \n",*R_activity_z); 
    fprintf(stderr,"  - N_attenuation_x:              %d \n",*N_attenuation_x);
    fprintf(stderr,"  - N_attenuation_y:              %d \n",*N_attenuation_y);
    fprintf(stderr,"  - N_attenuation_z:              %d \n",*N_attenuation_z);
    fprintf(stderr,"  - attenuation_size_x:           %f \n",*attenuation_size_x);
    fprintf(stderr,"  - attenuation_size_y:           %f \n",*attenuation_size_y);
    fprintf(stderr,"  - attenuation_size_z:           %f \n",*attenuation_size_z);
    fprintf(stderr,"  - T_attenuation_x:              %f \n",*T_attenuation_x); 
    fprintf(stderr,"  - T_attenuation_y:              %f \n",*T_attenuation_y); 
    fprintf(stderr,"  - T_attenuation_z:              %f \n",*T_attenuation_z); 
    fprintf(stderr,"  - R_attenuation_x:              %f \n",*R_attenuation_x); 
    fprintf(stderr,"  - R_attenuation_y:              %f \n",*R_attenuation_y); 
    fprintf(stderr,"  - R_attenuation_z:              %f \n",*R_attenuation_z); 
    fprintf(stderr,"  - N_axial:                      %d \n",*N_axial);
    fprintf(stderr,"  - N_azimuthal:                  %d \n",*N_azimuthal);  
    fprintf(stderr,"  - N_u:                          %d \n",*N_u);
    fprintf(stderr,"  - N_v:                          %d \n",*N_v);
    fprintf(stderr,"  - size_u:                       %f \n",*size_u);
    fprintf(stderr,"  - size_v:                       %f \n",*size_v);
    fprintf(stderr,"  - N_locations:                  %d \n",*N_locations);
    fprintf(stderr,"  - N_samples:                    %d \n",*N_samples); 
    fprintf(stderr,"  - sample_step:                  %f \n",*sample_step);
    fprintf(stderr,"  - background_activity:          %f \n",*background);
    fprintf(stderr,"  - background_attenuation:       %f \n",*background_attenuation); 
    fprintf(stderr,"  - truncate_negative_values:     %d \n",*truncate_negative_values);
    fprintf(stderr,"  - use_gpu:                      %d \n",*use_gpu); 
    fprintf(stderr,"  - offsets:                      %d %d %d %d %d %d %d %d  \n", offsets[0],offsets[1],offsets[2],offsets[3],offsets[4],offsets[5],offsets[6],offsets[7]); 
    fprintf(stderr,"  - locations:                    %d %d %d   %d %d %d   %d %d %d   %d %d %d   %d %d %d \n", locations[0],locations[1],locations[2],locations[3],locations[4],locations[5],locations[6],locations[7],locations[8],locations[9],locations[10],locations[11],locations[12],locations[13],locations[14]); 
    
}
    int status = STATUS_SUCCESS;
    int no_attenuation = 0;

    //attenuation or not?
    if (*N_attenuation_x == 0 && *N_attenuation_y == 0 && *N_attenuation_z == 0)
        no_attenuation = 1;

	// Allocate nifti images 
    int dim[8]; 

    // Allocate activity nifti images 
    dim[0]    = 3; 
    dim[1]    = *N_activity_x; 
    dim[2]    = *N_activity_y; 
    dim[3]    = *N_activity_z; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *activityImage = NULL; 
    activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    activityImage->data = (float *)(activity); 
    
    //  1) define sform matrix: 
    mat44 *scale_m       = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *rotation_m    = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *m             = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_scale_matrix(scale_m, (*activity_size_x)/(*N_activity_x), (*activity_size_y)/(*N_activity_y), (*activity_size_z)/(*N_activity_z));  
    et_create_translation_matrix(translation_m, -(*T_activity_x), -(*T_activity_y), -(*T_activity_z)); 
    et_create_rotation_matrix(rotation_m, *R_activity_x, *R_activity_y, *R_activity_z, 0, 0 ,0 , XYZ_ROTATION);
    *m = reg_mat44_mul(translation_m, scale_m);
	*m = reg_mat44_mul(rotation_m, m);

    //  2) set sform matrix: 
    activityImage->pixdim[0] = 1;    
    activityImage->qform_code = 0; 
    activityImage->sform_code = 1; 
    activityImage->sto_xyz.m[0][0]=m->m[0][0];     activityImage->sto_xyz.m[0][1]=m->m[0][1];    activityImage->sto_xyz.m[0][2]=m->m[0][2];     activityImage->sto_xyz.m[0][3]=m->m[0][3];
    activityImage->sto_xyz.m[1][0]=m->m[1][0];     activityImage->sto_xyz.m[1][1]=m->m[1][1];    activityImage->sto_xyz.m[1][2]=m->m[1][2];     activityImage->sto_xyz.m[1][3]=m->m[1][3];
    activityImage->sto_xyz.m[2][0]=m->m[2][0];     activityImage->sto_xyz.m[2][1]=m->m[2][1];    activityImage->sto_xyz.m[2][2]=m->m[2][2];     activityImage->sto_xyz.m[2][3]=m->m[2][3]; 
    activityImage->sto_xyz.m[3][0]=m->m[3][0];     activityImage->sto_xyz.m[3][1]=m->m[3][1];    activityImage->sto_xyz.m[3][2]=m->m[3][2];     activityImage->sto_xyz.m[3][3]=m->m[3][3];
    activityImage->sto_ijk = nifti_mat44_inverse(activityImage->sto_xyz); 
            
    // Allocate attenuation nifti image 
    nifti_image *attenuationImage = NULL;
    if(!no_attenuation)
            {
            dim[1]    = *N_attenuation_x; 
            dim[2]    = *N_attenuation_y; 
            dim[3]    = *N_attenuation_z; 
            dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);    
            // FIXME: set sform 
            }

    
    //Do projection
#ifdef _USE_CUDA
    if (*use_gpu)
        status = PET_project_compressed_gpu(activityImage, attenuationImage, projection, 
                                            offsets, locations, active, *N_locations, *N_axial, *N_azimuthal, 
                                            angles, *N_u, *N_v, *size_u, *size_v, 
                                            *N_samples, *sample_step, *background, *background_attenuation, *truncate_negative_values,
                                            *direction, *block_size, time_profiling);
    else
        status = PET_project_compressed_cpu(activityImage, attenuationImage, projection, 
                                        offsets, locations, active, *N_locations, *N_axial, *N_azimuthal, 
                                        angles, *N_u, *N_v, *size_u, *size_v, 
                                        *N_samples, *sample_step, *background, *background_attenuation, *truncate_negative_values,
                                        *direction); 
#else
    if (*use_gpu)
        fprintf_verbose( "PET_array_project_compressed: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
    status = PET_project_compressed_cpu(activityImage, attenuationImage, projection, 
                                    offsets, locations, active, *N_locations, *N_axial, *N_azimuthal, 
                                    angles, *N_u, *N_v, *size_u, *size_v, 
                                    *N_samples, *sample_step, *background, *background_attenuation, *truncate_negative_values,
                                    *direction); 
#endif

    //Free
    if( activityImage->fname != NULL ) free(activityImage->fname) ;
    if( activityImage->iname != NULL ) free(activityImage->iname) ;
    (void)nifti_free_extensions( activityImage ) ;
    free(activityImage) ;

    if(!no_attenuation)
        {
        if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
        if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
        (void)nifti_free_extensions( attenuationImage ) ;
        free(attenuationImage) ;
        }
    
    free(scale_m); 
    free(translation_m);
    free(rotation_m);
    free(m); 
    
    return status;	 
}




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
                           unsigned int *active )
{
    int status = STATUS_SUCCESS;
    int no_attenuation = 0;

    //attenuation or not?
    if (*N_attenuation_x == 0 && *N_attenuation_y == 0 && *N_attenuation_z == 0)
        no_attenuation = 1;
    
	// Allocate nifti images
    int dim[8]; 

    // Allocate activity nifti images
    dim[0]    = 3; 
    dim[1]    = *N_activity_x; 
    dim[2]    = *N_activity_y; 
    dim[3]    = *N_activity_z; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *activityImage = NULL;
    activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    activityImage->data = (float *)(activity);
    
    // Allocate attenuation nifti image 
    nifti_image *attenuationImage = NULL;
    if(!no_attenuation)
            {
            dim[0]    = 3; 
            dim[1]    = *N_attenuation_x; 
            dim[2]    = *N_attenuation_y; 
            dim[3]    = *N_attenuation_z; 
            dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);     
            }

    // Do projection: 
#ifdef _USE_CUDA
    status = PET_project_compressed_gpu_test(activityImage, attenuationImage, projection, 
               offsets, locations, active, *N_locations, *N_axial, *N_azimuthal); 
#endif 
    return status; 

    //Free:
    if( activityImage->fname != NULL ) free(activityImage->fname) ;
    if( activityImage->iname != NULL ) free(activityImage->iname) ;
    (void)nifti_free_extensions( activityImage ) ;
    free(activityImage) ;

    if(!no_attenuation)
        {
        if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
        if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
        (void)nifti_free_extensions( attenuationImage ) ;
        free(attenuationImage) ;
        }
    
    return status;	
}
                           

               
               





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
                           unsigned int *block_size, float *time_profiling )
{
    fprintf(stderr,"Compressed backprojection ..\n");
if (0) {
    fprintf(stderr,"  - N_activity_x:                 %d \n",*N_activity_x);
    fprintf(stderr,"  - N_activity_y:                 %d \n",*N_activity_y);
    fprintf(stderr,"  - N_activity_z:                 %d \n",*N_activity_z);
    fprintf(stderr,"  - activity_size_x:              %f \n",*activity_size_x);
    fprintf(stderr,"  - activity_size_y:              %f \n",*activity_size_y);
    fprintf(stderr,"  - activity_size_z:              %f \n",*activity_size_z);
    fprintf(stderr,"  - T_activity_x:                 %f \n",*T_activity_x); 
    fprintf(stderr,"  - T_activity_y:                 %f \n",*T_activity_y); 
    fprintf(stderr,"  - T_activity_z:                 %f \n",*T_activity_z); 
    fprintf(stderr,"  - R_activity_x:                 %f \n",*R_activity_x); 
    fprintf(stderr,"  - R_activity_y:                 %f \n",*R_activity_y); 
    fprintf(stderr,"  - R_activity_z:                 %f \n",*R_activity_z); 
    fprintf(stderr,"  - N_attenuation_x:              %d \n",*N_attenuation_x);
    fprintf(stderr,"  - N_attenuation_y:              %d \n",*N_attenuation_y);
    fprintf(stderr,"  - N_attenuation_z:              %d \n",*N_attenuation_z);
    fprintf(stderr,"  - attenuation_size_x:           %f \n",*attenuation_size_x);
    fprintf(stderr,"  - attenuation_size_y:           %f \n",*attenuation_size_y);
    fprintf(stderr,"  - attenuation_size_z:           %f \n",*attenuation_size_z);
    fprintf(stderr,"  - T_attenuation_x:              %f \n",*T_attenuation_x); 
    fprintf(stderr,"  - T_attenuation_y:              %f \n",*T_attenuation_y); 
    fprintf(stderr,"  - T_attenuation_z:              %f \n",*T_attenuation_z); 
    fprintf(stderr,"  - R_attenuation_x:              %f \n",*R_attenuation_x); 
    fprintf(stderr,"  - R_attenuation_y:              %f \n",*R_attenuation_y); 
    fprintf(stderr,"  - R_attenuation_z:              %f \n",*R_attenuation_z); 
    fprintf(stderr,"  - N_axial:                      %d \n",*N_axial);
    fprintf(stderr,"  - N_azimuthal:                  %d \n",*N_azimuthal); 
    fprintf(stderr,"  - N_u:                          %d \n",*N_u);
    fprintf(stderr,"  - N_v:                          %d \n",*N_v);
    fprintf(stderr,"  - size_u:                       %f \n",*size_u);
    fprintf(stderr,"  - size_v:                       %f \n",*size_v);
    fprintf(stderr,"  - N_locations:                  %d \n",*N_locations);
    fprintf(stderr,"  - use_gpu:                      %d \n",*use_gpu); 
    fprintf(stderr,"  - N_samples:                    %d \n",*N_samples); 
    fprintf(stderr,"  - sample_step:                  %f \n",*sample_step); 
    fprintf(stderr,"  - background_activity:          %f \n",*background_activity); 
    fprintf(stderr,"  - background_attenuation:       %f \n",*background_attenuation); 
    fprintf(stderr,"  - direction:                    %d \n",*direction); 
    fprintf(stderr,"  - block_size:                   %d \n",*block_size); 
}
               
    int status = STATUS_SUCCESS;
    int no_attenuation = 0;

    //attenuation or not?
    if (*N_attenuation_x == 0 && *N_attenuation_y == 0 && *N_attenuation_z == 0)
        no_attenuation = 1;

	// Allocate nifti images 
    int dim[8]; 

    // Allocate back-projection nifti image
    dim[0]    = 3; 
    dim[1]    = *N_activity_x; 
    dim[2]    = *N_activity_y; 
    dim[3]    = *N_activity_z; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *backprojectionImage = NULL; 
    backprojectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    backprojectionImage->data = (float *)(backprojection); 
    
    //  1) define sform matrix: 
    mat44 *scale_m       = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *rotation_m    = (mat44 *)calloc(1,sizeof(mat44)); 
    mat44 *m             = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_scale_matrix(scale_m, (*activity_size_x)/(*N_activity_x), (*activity_size_y)/(*N_activity_y), (*activity_size_z)/(*N_activity_z));  
    et_create_translation_matrix(translation_m, -(*T_activity_x), -(*T_activity_y), -(*T_activity_z)); 
    et_create_rotation_matrix(rotation_m, *R_activity_x, *R_activity_y, *R_activity_z, 0, 0 ,0 , XYZ_ROTATION);
    *m = reg_mat44_mul(translation_m, scale_m);
	*m = reg_mat44_mul(rotation_m, m);

    //  2) set sform matrix: 
    backprojectionImage->pixdim[0] = 1;    
    backprojectionImage->qform_code = 0; 
    backprojectionImage->sform_code = 1; 
    backprojectionImage->sto_xyz.m[0][0]=m->m[0][0];     backprojectionImage->sto_xyz.m[0][1]=m->m[0][1];    backprojectionImage->sto_xyz.m[0][2]=m->m[0][2];     backprojectionImage->sto_xyz.m[0][3]=m->m[0][3];
    backprojectionImage->sto_xyz.m[1][0]=m->m[1][0];     backprojectionImage->sto_xyz.m[1][1]=m->m[1][1];    backprojectionImage->sto_xyz.m[1][2]=m->m[1][2];     backprojectionImage->sto_xyz.m[1][3]=m->m[1][3];
    backprojectionImage->sto_xyz.m[2][0]=m->m[2][0];     backprojectionImage->sto_xyz.m[2][1]=m->m[2][1];    backprojectionImage->sto_xyz.m[2][2]=m->m[2][2];     backprojectionImage->sto_xyz.m[2][3]=m->m[2][3]; 
    backprojectionImage->sto_xyz.m[3][0]=m->m[3][0];     backprojectionImage->sto_xyz.m[3][1]=m->m[3][1];    backprojectionImage->sto_xyz.m[3][2]=m->m[3][2];     backprojectionImage->sto_xyz.m[3][3]=m->m[3][3];
    backprojectionImage->sto_ijk = nifti_mat44_inverse(backprojectionImage->sto_xyz); 
            
    // Allocate attenuation nifti image 
    nifti_image *attenuationImage = NULL;
    if(!no_attenuation)
            {
            dim[1]    = *N_attenuation_x; 
            dim[2]    = *N_attenuation_y; 
            dim[3]    = *N_attenuation_z; 
            dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);    
            // FIXME: set sform 
            }

    
    //Do projection
#ifdef _USE_CUDA
    if (*use_gpu)
        status = PET_backproject_compressed_gpu(backprojectionImage, attenuationImage, projection_data, 
                                            offsets, locations, active, *N_locations, *N_axial, *N_azimuthal, 
                                            angles, *N_u, *N_v, *size_u, *size_v, 
                                            *N_samples, *sample_step, *background_activity, *background_attenuation, 
                                            *direction, *block_size, time_profiling);
    else
        status = PET_backproject_compressed_cpu(backprojectionImage, attenuationImage, projection_data, 
                                        offsets, locations, active, *N_locations, *N_axial, *N_azimuthal, 
                                        angles, *N_u, *N_v, *size_u, *size_v, 
                                        *N_samples, *sample_step, *background_activity, *background_attenuation, 
                                        *direction); 
#else
    if (*use_gpu)
        fprintf_verbose( "PET_array_project_compressed: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
    status = PET_backproject_compressed_cpu(backprojectionImage, attenuationImage, projection_data, 
                                    offsets, locations, active, *N_locations, *N_axial, *N_azimuthal, 
                                    angles, *N_u, *N_v, *size_u, *size_v, 
                                    *N_samples, *sample_step, *background_activity, *background_attenuation, 
                                    *direction); 
#endif

    //Free
    if( backprojectionImage->fname != NULL ) free(backprojectionImage->fname) ;
    if( backprojectionImage->iname != NULL ) free(backprojectionImage->iname) ;
    (void)nifti_free_extensions( backprojectionImage ) ;
    free(backprojectionImage) ;

    if(!no_attenuation)
        {
        if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
        if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
        (void)nifti_free_extensions( attenuationImage ) ;
        free(attenuationImage) ;
        }
    
    free(scale_m); 
    free(translation_m);
    free(rotation_m);
    free(m); 
    
    return status;	 
}







nifti_image *ET_create_nifti_3D_scalar(unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *pixdim_x, float *pixdim_y, float *pixdim_z, float *image_array, float *affine) 
{
    int dim[8]; 
    dim[0]    = 3; 
    dim[1]    = *Nx; 
    dim[2]    = *Ny; 
    dim[3]    = *Nz; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *nii = NULL;
    nii = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    nii->data = (float *) image_array;
    nii->pixdim[0] = 1;
    nii->pixdim[1] = *pixdim_x; 
    nii->pixdim[2] = *pixdim_y; 
    nii->pixdim[3] = *pixdim_z; 
    nii->dx = *pixdim_x; 
    nii->dy = *pixdim_y;
    nii->dz = *pixdim_z;
    // set the sform affine matrix 
    if ( affine != NULL) { 
        nii->qform_code = 0; 
        nii->sform_code = 1; 
        nii->sto_xyz.m[0][0]=affine[0];     nii->sto_xyz.m[0][1]=affine[1];    nii->sto_xyz.m[0][2]=affine[2];     nii->sto_xyz.m[0][3]=affine[3];
        nii->sto_xyz.m[1][0]=affine[4];     nii->sto_xyz.m[1][1]=affine[5];    nii->sto_xyz.m[1][2]=affine[6];     nii->sto_xyz.m[1][3]=affine[7];
        nii->sto_xyz.m[2][0]=affine[8];     nii->sto_xyz.m[2][1]=affine[9];    nii->sto_xyz.m[2][2]=affine[10];    nii->sto_xyz.m[2][3]=affine[11]; 
        nii->sto_xyz.m[3][0]=affine[12];    nii->sto_xyz.m[3][1]=affine[13];   nii->sto_xyz.m[3][2]=affine[14];    nii->sto_xyz.m[3][3]=affine[15];
        nii->sto_ijk = nifti_mat44_inverse(nii->sto_xyz); 
        }
    else {
        nii->qform_code = 0; 
        nii->sform_code = 0;   
    }
    return nii; 
}

nifti_image *ET_create_nifti_3D_vector(unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *dim_vector, float *pixdim_x, float *pixdim_y, float *pixdim_z, float *image_array, float *affine) 
{
    int dim[8]; 
    dim[0]    = 5; 
    dim[1]    = *Nx; 
    dim[2]    = *Ny; 
    dim[3]    = *Nz; 
    dim[4]    = 1; 
    dim[5]    = *dim_vector; 
    dim[6]=1; dim[7]=1; 
    nifti_image *nii = NULL;
    nii = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    nii->data = (float *) image_array;
    nii->pixdim[0] = 1;
    nii->pixdim[1] = *pixdim_x; 
    nii->pixdim[2] = *pixdim_y; 
    nii->pixdim[3] = *pixdim_z; 
    nii->pixdim[4] = 1.0; nii->pixdim[5] = 1.0; nii->pixdim[6] = 1.0; nii->pixdim[7] = 1.0;
    nii->dx=*pixdim_x; nii->dy=*pixdim_y; nii->dz=*pixdim_z; 
    nii->dt=1.0; nii->du=1.0; nii->dv=1.0; nii->dw=1.0; 
    nii->nx=*Nx; nii->ny=*Ny; nii->nz=*Nz; 
    nii->nt=1;   nii->nu=3;   nii->nv=1;   nii->nw=1; 
    nii->nvox=nii->nx*nii->ny*nii->nz*nii->nt; //*nii->nu; 
    // set the sform affine matrix 
    if ( affine != NULL) { 
        nii->qform_code = 0; 
        nii->sform_code = 1; 
        nii->sto_xyz.m[0][0]=affine[0];     nii->sto_xyz.m[0][1]=affine[1];    nii->sto_xyz.m[0][2]=affine[2];     nii->sto_xyz.m[0][3]=affine[3];
        nii->sto_xyz.m[1][0]=affine[4];     nii->sto_xyz.m[1][1]=affine[5];    nii->sto_xyz.m[1][2]=affine[6];     nii->sto_xyz.m[1][3]=affine[7];
        nii->sto_xyz.m[2][0]=affine[8];     nii->sto_xyz.m[2][1]=affine[9];    nii->sto_xyz.m[2][2]=affine[10];    nii->sto_xyz.m[2][3]=affine[11]; 
        nii->sto_xyz.m[3][0]=affine[12];    nii->sto_xyz.m[3][1]=affine[13];   nii->sto_xyz.m[3][2]=affine[14];    nii->sto_xyz.m[3][3]=affine[15];
        nii->sto_ijk = nifti_mat44_inverse(nii->sto_xyz); 
        }
    else {
        nii->qform_code = 0; 
        nii->sform_code = 0;   
    }
    return nii; 
}


#define FREE_DATA    1
#define NO_FREE_DATA 0

void ET_free_nifti(nifti_image* nii, unsigned int free_data)
{
    if( free_data==FREE_DATA ) if (nii->data != NULL) free(nii->data); 
    if( nii->fname != NULL ) free(nii->fname) ;
    if( nii->iname != NULL ) free(nii->iname) ;
    (void)nifti_free_extensions( nii ) ;
    free(nii) ;
}


unsigned int ET_identity_transform_44(float *array)
{
    for (int i=0; i<16; i++)
         array[i] = 0; 
    array[0] = 1; array[5]=1; array[10]=1; array[15]=1; 
    return STATUS_SUCCESS; 
}

unsigned int ET_combine_status(unsigned int new_status, unsigned int old_status)
{
    // FIXME: do or of new_status and old_status. Define status flags accordingly: one bit per type, 0 for success. 
    return new_status; 
}


unsigned int ET_spherical_phantom(float *image_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *sizex, float *sizey, float *sizez, float *centerx, float *centery, float *centerz, float *radius, float *inner_value, float *outer_value)
{
    unsigned int status=STATUS_SUCCESS; 
    // Allocate the nifti image
    float pixdim_x = (*sizex)/(*Nx); float pixdim_y = (*sizex)/(*Nx); float pixdim_z = (*sizex)/(*Nx);
    nifti_image *phantomImage = ET_create_nifti_3D_scalar(Nx, Ny, Nz, &pixdim_x, &pixdim_y, &pixdim_z, image_array, NULL); 
    
    // Make spherical phantom 
    status = et_spherical_phantom(phantomImage, *centerx, *centery, *centerz, *radius, *inner_value, *outer_value); 

    //Free
    ET_free_nifti( phantomImage, NO_FREE_DATA); 

    return status; 
}


unsigned int ET_cylindrical_phantom(float *image_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *sizex, float *sizey, float *sizez, float *centerx, float *centery, float *centerz, float *radius, float *length, unsigned int *axis, float *inner_value, float *outer_value)
{
    unsigned int status=STATUS_SUCCESS; 
    // Allocate the nifti image
    float pixdim_x = (*sizex)/(*Nx); float pixdim_y = (*sizex)/(*Nx); float pixdim_z = (*sizex)/(*Nx);
    nifti_image *phantomImage = ET_create_nifti_3D_scalar(Nx, Ny, Nz, &pixdim_x, &pixdim_y, &pixdim_z, image_array, NULL); 
    
    // Make spherical phantom 
    status = et_cylindrical_phantom(phantomImage, *centerx, *centery, *centerz, *radius, *length, *axis, *inner_value, *outer_value); 

    //Free
    ET_free_nifti(phantomImage, NO_FREE_DATA ); 

    return status; 
}


unsigned int ET_spheres_ring_phantom(float *image_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *sizex, float *sizey, float *sizez, float *centerx, float *centery, float *centerz, float *ring_radius, float *min_sphere_radius, float *max_sphere_radius, unsigned int *N_spheres, float *inner_value, float *outer_value, float *taper, unsigned int *ring_axis)
{
    unsigned int status=STATUS_SUCCESS; 
    // Allocate the nifti image
    float pixdim_x = (*sizex)/(*Nx); float pixdim_y = (*sizex)/(*Nx); float pixdim_z = (*sizex)/(*Nx);
    nifti_image *phantomImage = ET_create_nifti_3D_scalar(Nx, Ny, Nz, &pixdim_x, &pixdim_y, &pixdim_z, image_array, NULL); 
    
    // Make phantom 
    status = et_spheres_ring_phantom(phantomImage, *centerx, *centery, *centerz, *ring_radius, *min_sphere_radius, *max_sphere_radius, *N_spheres, *inner_value, *outer_value, *taper, *ring_axis); 

    //nifti_set_filenames(phantomImage,"/Users/spedemon/Desktop/phantom.nii",0,0); 
    //nifti_image_write(phantomImage); 
    
    //Free
    ET_free_nifti( phantomImage, NO_FREE_DATA); 

    return status; 
}


unsigned int TR_resample_grid(float *resampled_image_array, float *image_array, float *affine, float *grid_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode)
{
    unsigned int global_status = STATUS_SUCCESS; 
    unsigned int running_status = STATUS_SUCCESS; 

    // allocate nifti image for input image and set affine transformation
    float pixdim_x = 1; float pixdim_y = 1; float pixdim_z = 1;
    nifti_image *image = ET_create_nifti_3D_scalar(Nx, Ny, Nz, &pixdim_x, &pixdim_y, &pixdim_z, image_array, affine); 
    
    // allocate nifti image for resample grid
    float *identity = (float *)malloc(16*sizeof(float)); 
    ET_identity_transform_44(identity); 
    unsigned int dim_vector = 3; 
    nifti_image *grid = ET_create_nifti_3D_vector(Nx_grid, Ny_grid, Nz_grid, &dim_vector, &pixdim_x, &pixdim_y, &pixdim_z, grid_array, identity); 
    
    // allocate nifti image for the resampled image 
    nifti_image *resampled = ET_create_nifti_3D_scalar(Nx_grid, Ny_grid, Nz_grid, &pixdim_x, &pixdim_y, &pixdim_z, resampled_image_array, identity );
    
    // resample the image 
	#ifdef _USE_CUDA
    if ( *use_gpu )
        running_status = tr_resample_grid_gpu(resampled, image, grid, *background, *interpolation_mode); 
    else
    #endif 
        running_status = tr_resample_grid_cpu(resampled, image, grid, *background, *interpolation_mode); 
    global_status = ET_combine_status(running_status, global_status); 
    
    // free memory
    free(identity); 
    ET_free_nifti( image, NO_FREE_DATA ); 
    ET_free_nifti( resampled, NO_FREE_DATA ); 
    ET_free_nifti( grid, NO_FREE_DATA ); 
     
    return global_status; 
}


unsigned int TR_resample_box(float *resampled_image_array, float *image_array, float *affine, float *min_coords, float *max_coords, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode)
{
    unsigned int status = STATUS_SUCCESS; 
    return status; 
}

unsigned int TR_gradient_grid(float *gradient_array, float *image_array, float *affine, float *grid, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode)
{
    unsigned int status = STATUS_SUCCESS; 
    return status; 
}

unsigned int TR_gradient_box(float *gradient_array, float *image_array, float *affine, float *min_coords, float *max_coords, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, unsigned int *Nx_grid, unsigned int *Ny_grid, unsigned int *Nz_grid, float *background, unsigned int *use_gpu, unsigned int *interpolation_mode)
{
    unsigned int status = STATUS_SUCCESS; 
    return status; 
}



unsigned int TR_grid_from_box_and_affine(float *grid_array, float *affine_box2grid, float *min_x, float *min_y, float *min_z, float *max_x, float *max_y, float *max_z, unsigned int *n_x, unsigned int * n_y, unsigned int *n_z)
{
    unsigned int global_status = STATUS_SUCCESS; 
    unsigned int running_status = STATUS_SUCCESS; 

    // instantiate nifti image for the grid and set affine (i.e. sform)

    unsigned int dim_vector = 3; 
    float        pixdim_x = ((*max_x) - (*min_x))/(*n_x - 1); 
    float        pixdim_y = ((*max_y) - (*min_y))/(*n_y - 1); 
    float        pixdim_z = ((*max_z) - (*min_z))/(*n_z - 1); 

    nifti_image *grid = ET_create_nifti_3D_vector(n_x, n_y, n_z, &dim_vector, &pixdim_x, &pixdim_y, &pixdim_z, grid_array, affine_box2grid); 
    
    mat44 *scale = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_scale_matrix(scale, grid->pixdim[1], grid->pixdim[2], grid->pixdim[3]); 
    mat44 *translation = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_translation_matrix(translation, *min_x, *min_y, *min_z); 
    mat44 m = reg_mat44_mul(translation, scale); 
    
    // compute the position of the grid points 

    running_status = ET_box_to_grid(grid, &m); 
    global_status = ET_combine_status(running_status, global_status); 

    // free memory 

    ET_free_nifti( grid, NO_FREE_DATA); 
    free(scale); 
    free(translation); 
    
    return global_status; 
}



unsigned int TR_transform_grid(float *transformed_grid_array, float *grid_array, unsigned int *Nx, unsigned int *Ny, unsigned int *Nz, float *affine_from_grid, unsigned int *use_gpu)
{
    unsigned int running_status = STATUS_SUCCESS; 
    unsigned int global_status  = STATUS_SUCCESS; 
    
    // 1) instantiate nifti_image for grid
    unsigned int dim_vector = 3; 
    float        pixdim_x = 1.0; 
    float        pixdim_y = 1.0; 
    float        pixdim_z = 1.0; 

    nifti_image *grid = ET_create_nifti_3D_vector(Nx, Ny, Nz, &dim_vector, &pixdim_x, &pixdim_y, &pixdim_z, grid_array, NULL); 

    // 2) instantiate nifti_image for transformed image 
    nifti_image *transformed_grid = ET_create_nifti_3D_vector(Nx, Ny, Nz, &dim_vector, &pixdim_x, &pixdim_y, &pixdim_z, transformed_grid_array, NULL); 

    // 3) transform 
    mat44 *affine = (mat44 *)calloc(1,sizeof(mat44)); 
    et_create_affine_matrix_from_buffer(affine, affine_from_grid); 

    #ifdef _USE_CUDA
    if (*use_gpu)
        running_status = tr_transform_grid_gpu(transformed_grid, grid, affine); 
    else
    #endif
        running_status = tr_transform_grid_cpu(transformed_grid, grid, affine); 
    global_status = ET_combine_status(running_status, global_status); 

    // free memory 
    free(affine); 
    ET_free_nifti( grid, NO_FREE_DATA); 
    ET_free_nifti( transformed_grid, NO_FREE_DATA); 
    
    return global_status; 
}





// Compress and uncompress projection data 
extern "C" unsigned int PET_initialize_compression_structure(unsigned int *N_axial, unsigned int *N_azimuthal, unsigned int *N_u, unsigned int *N_v, int *offsets_matrix, unsigned short *locations )
{
    return pet_initialize_compression_structure(N_axial, N_azimuthal, N_u, N_v, offsets_matrix, locations ); 
}

extern "C" unsigned int PET_compress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection)
{
    return pet_compress_projection(n_locations, n_axial, n_azimuthal, n_u, n_v, offsets_matrix, counts, locations, projection); 
}

extern "C" unsigned int PET_uncompress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection)
{
    return pet_uncompress_projection(n_locations, n_axial, n_azimuthal, n_u, n_v, offsets_matrix, counts, locations, projection); 
}

extern "C" unsigned int PET_compress_projection_array(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, float *threshold, unsigned int *N_active, bool *active, float *data, int *offsets, unsigned short *locations, float *compressed_data)
{
    return pet_compress_projection_array(N_ax, N_az, N_u, N_v, threshold, N_active, active, data, offsets, locations, compressed_data); 
}

extern "C" unsigned int PET_get_subset_sparsity(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, unsigned int *N_locations, int *offsets, unsigned short *locations, unsigned int *subsets_matrix, int *offsets_sub, unsigned short *locations_sub) 
{
    return pet_get_subset_sparsity(N_ax, N_az, N_u, N_v, N_locations, offsets, locations, subsets_matrix, offsets_sub, locations_sub);
}

extern "C" unsigned int PET_get_subset_projection_array(unsigned int *N_ax, unsigned int *N_az, unsigned int *N_u, unsigned int *N_v, unsigned int *N_locations, float *data, int *offsets, unsigned short *locations, int *offsets_sub, unsigned short *locations_sub, unsigned int *subsets_matrix, float *data_sub)
{
    return pet_get_subset_projection_array(N_ax, N_az, N_u, N_v, N_locations, data, offsets, locations, offsets_sub, locations_sub, subsets_matrix, data_sub); 
}



