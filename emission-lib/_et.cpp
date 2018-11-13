/*
 *  _et.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Harvard University, Martinos Center for Biomedical Imaging
 *  Jan. 2014.
 */


#include "_et.h"
#include "_et_common.h"
#include <stdio.h>
#include <sys/time.h>
#include <stdbool.h>

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define PI 3.141592653589793


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   CPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned short int et_is_block_multiple(unsigned short int size)
{
    if (size % (ET_BLOCK_SIZE*ET_BLOCK_SIZE) == 0)
        return 1;
    return 0;
}
unsigned short int et_get_block_size(void)
{
    return ET_BLOCK_SIZE*ET_BLOCK_SIZE;
}


//! Affine transformation of nifti_image 
/*!
  \param *sourceImage the source image to be transformed. 
  \param *transformedImage the transformed image. 
  \param *affineTransformation the [4x4] transformed matrix. 
  \param background the background value when resampling the transformed image. 
*/
int et_affine(nifti_image *sourceImage, nifti_image *transformedImage, mat44 *affineTransformation, float background)
{
    /* Allocate the deformation Field image */
    nifti_image *positionFieldImage = nifti_copy_nim_info(sourceImage);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[1]=positionFieldImage->nx=sourceImage->nx;
    positionFieldImage->dim[2]=positionFieldImage->ny=sourceImage->ny;
    positionFieldImage->dim[3]=positionFieldImage->nz=sourceImage->nz;
    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
    positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    positionFieldImage->nbyper = sizeof(float);
    positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
    
    /* Apply affine */
    reg_affine_positionField(       affineTransformation,
                                    sourceImage,
                                    positionFieldImage );
    /* Resample the source image */
    reg_resampleSourceImage<float>( sourceImage,
                                    sourceImage,
                                    transformedImage,
                                    positionFieldImage,
                                    NULL,
                                    1,
                                    background);
    nifti_image_free(positionFieldImage);
    return 0;
}


//! Rotate a nifti_image in 3D 
/*!
  \param *sourceImage the source image to be transformed. 
  \param *resultImage the rotated image. 
  \param theta_x the rotation angle around x axis in radians. 
  \param theta_y the rotation angle around y axis in radians. 
  \param theta_z the rotation angle around z axis in radians. 
  \param center_x the center of rotation along x.  
  \param center_y the center of rotation along y. 
  \param center_z the center of rotation along z. 
  \param background the background value when resampling the transformed image. 
*/
int et_rotate(nifti_image *sourceImage, nifti_image *resultImage, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, float background, int axis_order)
{
    int status; 
    //Create transformation matrix
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    et_create_rotation_matrix(affineTransformation, theta_x, theta_y, theta_z, center_x, center_y, center_z, axis_order);
    
    //Apply affine transformation
    status = et_affine(sourceImage, resultImage, affineTransformation, background);

    //Free
    free(affineTransformation);

    return status;
}



//! Projection for Emission Imaging
/*!
  \param *activityImage the activity (or its estimate). NULL for attenuation and background activity only. 
  \param *sinoImage the photon counts in projection space. 
  \param *psfImage the depth-dependent point spread function, NULL for no PSF. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_project(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
        int separable_psf = 0;
        int psf_size[3];

        /* Check consistency of input */
        nifti_image *referenceImage;              // this image holds information about image size and voxel size (activity or attenuation might not be defined (NULL pointers) ) 
        if (activityImage==NULL && attenuationImage==NULL)
            {
            fprintf(stderr, "et_project: Error - define at least one between activityImage and attenuationImage. \n");
            return niftyrec_error_parameters; 
            }
        else if (attenuationImage==NULL)
            referenceImage=activityImage;
        else
            referenceImage=attenuationImage; 

    /* Allocate the deformation Field image */
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);
     
    nifti_image *positionFieldImage = nifti_copy_nim_info(referenceImage);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[1]=positionFieldImage->nx=referenceImage->nx;
    positionFieldImage->dim[2]=positionFieldImage->ny=referenceImage->ny;
    positionFieldImage->dim[3]=positionFieldImage->nz=referenceImage->nz;
    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
    positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    positionFieldImage->nbyper = sizeof(float);
    positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
        if (positionFieldImage->data==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
        alloc_record_add(memory_record,(void*)positionFieldImage,ALLOCTYPE_NIFTI);
    
    /* Allocate arrays */
        int dim[8];
    dim[0]    = 3;
    dim[1]    = referenceImage->dim[1];
    dim[2]    = referenceImage->dim[2];
    dim[3]    = referenceImage->dim[3];
    dim[4]    = 1;
    dim[5]    = 1;
    dim[6]    = 1;
    dim[7]    = 1;
    nifti_image *rotatedImage=NULL;
        if (activityImage!=NULL)
            {
            rotatedImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedImage->data = (float *)malloc(referenceImage->nvox*sizeof(float));    
            if (rotatedImage->data==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
            alloc_record_add(memory_record,(void*)rotatedImage,ALLOCTYPE_NIFTI);
            }

    nifti_image *rotatedAttenuationImage=NULL;
        if (attenuationImage != NULL)
            {
            rotatedAttenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedAttenuationImage->data = (float *)malloc(referenceImage->nvox*sizeof(float));
            if (rotatedAttenuationImage->data==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
            alloc_record_add(memory_record,(void*)rotatedAttenuationImage,ALLOCTYPE_NIFTI);    
            }
        
    /* Alloc transformation matrix */
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        if (affineTransformation==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
        alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_HOST);

        /* Decide whether to use FFT or separate the convolution */
        float *psfSeparated=NULL;
        float psf_norm;
    if(psfImage!=NULL)
            {
            if(1) //(psfImage->nx <= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1) 
                {
                separable_psf=1;
                psf_size[0] = psfImage->dim[1];
                psf_size[1] = psfImage->dim[2];
                psf_size[2] = psfImage->dim[3];
                psfSeparated = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
                if (psfSeparated==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
                alloc_record_add(memory_record,(void*)psfSeparated,ALLOCTYPE_HOST);
                for (int n=0; n<psf_size[2];n++) 
                    {
                    psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                    for (int i=0;i<psf_size[0];i++) 
                        {
                        psfSeparated[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                        psfSeparated[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                        }
                    }
                }
            }

        /* Project */
    for(int cam=0; cam<n_cameras; cam++){
        // Apply affine //
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
                fprintf_verbose( "et_project: Center:   %f  %f  %f  \n",centers[0*n_cameras+cam], centers[1*n_cameras+cam], centers[2*n_cameras+cam]);
        et_create_rotation_matrix(      affineTransformation, 
                        cameras[0*n_cameras+cam], 
                        cameras[1*n_cameras+cam], 
                        cameras[2*n_cameras+cam], 
                        centers[0*n_cameras+cam], 
                        centers[1*n_cameras+cam], 
                        centers[2*n_cameras+cam], 
                        XYZ_ROTATION);

        reg_affine_positionField(    affineTransformation,
                        referenceImage,
                        positionFieldImage );
        // Resample the source image //
                if (activityImage != NULL)
                    reg_resampleSourceImage<float>(activityImage,
                        activityImage,
                        rotatedImage,
                        positionFieldImage,
                        NULL,
                        1,
                        background
                        );    

                // Resample the attenuation map //
                if (attenuationImage != NULL)
                    {
                    reg_resampleSourceImage<float>(attenuationImage,
                        attenuationImage,
                        rotatedAttenuationImage,
                        positionFieldImage,
                        NULL,
                        1,
                        background_attenuation
                        );                        
                    }

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL && activityImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D( rotatedImage,
                                                psfSeparated,
                                                psf_size[0],
                                                psf_size[1],
                                                rotatedImage, 
                                                background );
                    else
                        et_convolve2D(          rotatedImage,
                                                psfImage,
                                                rotatedImage, 
                                                background );
                    }

        // Integrate along lines //
                et_line_integral_attenuated(    rotatedImage, 
                                                rotatedAttenuationImage, 
                                                sinoImage, 
                                                cam, 
                                                background); 
    }

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<sinoImage->nvox; i++)
                {
                if (sino_data[i] < 0)
                    sino_data[i] = 0;
                }
            }

    /* Deallocate memory */
        return alloc_record_destroy(memory_record); 
}



//! Back-projection for Emission Imaging
/*!
  \param *sinogramImage the data to be back-projected in projection space. 
  \param *backprojectionImage the output backprojection. 
  \param *psfImage the depth-dependent point spread function, NULL for no point spread function. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_backproject(nifti_image *sinogramImage, nifti_image *backprojectionImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{

        int separable_psf = 0;
        int psf_size[3];

    /* Allocate the deformation Field image */
        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);

    nifti_image *positionFieldImage = nifti_copy_nim_info(backprojectionImage);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[1]=positionFieldImage->nx=backprojectionImage->nx;
    positionFieldImage->dim[2]=positionFieldImage->ny=backprojectionImage->ny;
    positionFieldImage->dim[3]=positionFieldImage->nz=backprojectionImage->nz;
    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
    positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    positionFieldImage->nbyper = sizeof(float); 
    positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper); 
        if (positionFieldImage->data==NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,positionFieldImage,ALLOCTYPE_NIFTI);
    
    /* Allocate arrays */
        int dim[8];
    dim[0]    = 3;
    dim[1]    = backprojectionImage->dim[1];
    dim[2]    = backprojectionImage->dim[2];
    dim[3]    = backprojectionImage->dim[3];
    dim[4]    = 1;
    dim[5]    = 1;
    dim[6]    = 1;
    dim[7]    = 1;
    nifti_image *rotatedImage;
        rotatedImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        rotatedImage->data = (int *)malloc(backprojectionImage->nvox*sizeof(int));    
        if (rotatedImage->data == NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,rotatedImage,ALLOCTYPE_NIFTI); 

    nifti_image *temp_backprojectionImage;
        temp_backprojectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        temp_backprojectionImage->data = (int *)malloc(backprojectionImage->nvox*sizeof(int));            
        if (temp_backprojectionImage->data == NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,temp_backprojectionImage,ALLOCTYPE_NIFTI); 

    nifti_image *rotatedAttenuationImage;
        if (attenuationImage != NULL)
            {
            rotatedAttenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedAttenuationImage->data = (int *)malloc(attenuationImage->nvox*sizeof(int));    
            if (rotatedAttenuationImage->data == NULL) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_alloccpu; }
            alloc_record_add(memory_record,rotatedAttenuationImage,ALLOCTYPE_NIFTI); 
            }

    /* Clear accumulator */
    et_clear_accumulator(backprojectionImage);    
    
    /* Alloc transformation matrix */
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        if (affineTransformation == NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,affineTransformation,ALLOCTYPE_HOST); 

        /* Decide whether to use FFT or separate the convolution */
        float *psfSeparated=NULL;
        float psf_norm;
    if(psfImage!=NULL)
            {
            if(1) //(psfImage->nx <= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1) 
                {
                separable_psf=1;
                psf_size[0] = psfImage->dim[1];
                psf_size[1] = psfImage->dim[2];
                psf_size[2] = psfImage->dim[3];
                psfSeparated = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
                if (psfSeparated == NULL) {
                    alloc_record_destroy(memory_record); 
                    return niftyrec_error_alloccpu; }
                alloc_record_add(memory_record,psfSeparated,ALLOCTYPE_HOST); 
                for (int n=0; n<psf_size[2];n++) 
                    {
                    psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                    for (int i=0;i<psf_size[0];i++) 
                        {
                        psfSeparated[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                        psfSeparated[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                        }
                    }
                }
            }       

        /* Backproject */
    for(int cam=0; cam<n_cameras; cam++){
                /* Rotate attenuation */

                if (attenuationImage != NULL)                
                    {
                    et_create_rotation_matrix(    affineTransformation,
                        cameras[0*n_cameras+cam],
                        cameras[1*n_cameras+cam],
                        cameras[2*n_cameras+cam],
                        centers[0*n_cameras+cam],
                        centers[1*n_cameras+cam],
                            centers[2*n_cameras+cam],
                        XYZ_ROTATION);
                    reg_affine_positionField(    affineTransformation,
                        attenuationImage,
                        positionFieldImage);
                    reg_resampleSourceImage(    attenuationImage,
                        attenuationImage,
                        rotatedAttenuationImage,
                        positionFieldImage,
                        NULL,
                        1,
                        background_attenuation
                        );
                    }

        /* Line Backproject */
                if (attenuationImage != NULL)
                    {
                    et_line_backproject_attenuated(    sinogramImage,
                                                temp_backprojectionImage, 
                                                rotatedAttenuationImage, 
                                                cam );
                    }
                else
                    {
                    et_line_backproject(        sinogramImage,
                                                temp_backprojectionImage,
                                                cam );
                    }

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D( temp_backprojectionImage,
                                                psfSeparated,
                                                psf_size[0],
                                                psf_size[1],
                                                temp_backprojectionImage, 
                                                0.0f );
                    else
                        et_convolve2D(          temp_backprojectionImage,
                                                psfImage,
                                                temp_backprojectionImage, 
                                                0.0f );
                    }

        /* Rotate backprojection */
        et_create_rotation_matrix(    affineTransformation,
                        -cameras[0*n_cameras+cam],
                        -cameras[1*n_cameras+cam],
                        -cameras[2*n_cameras+cam],
                        centers[0*n_cameras+cam],
                        centers[1*n_cameras+cam],
                            centers[2*n_cameras+cam],
                        ZYX_ROTATION);
                        
        reg_affine_positionField(    affineTransformation,
                        backprojectionImage,
                        positionFieldImage);

        reg_resampleSourceImage(    temp_backprojectionImage,
                        temp_backprojectionImage,
                        rotatedImage,
                        positionFieldImage,
                        NULL,
                        1,
                        background
                        );

        /* Accumulate */
        et_accumulate(            rotatedImage,
                        backprojectionImage );
    }

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* accumulator_data = (float*) backprojectionImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<backprojectionImage->nvox; i++)
                {
                if (accumulator_data[i] < 0)
                   accumulator_data[i] = 0;
                }
            }
    /*Free*/
        return alloc_record_destroy(memory_record); 
}



//! Fisher Information Matrix of a grid of voxels, Emission Imaging
/*!
  \param from_projection whether the input image is a projection image (1) or activity image (0)
  \param *inputImage input image: projection image or activity image. 
  \param *gridImage grid image (same size as the activity), indexes from 1 to N_points at grid locations, 0 elsewhere. 
  \param *fisherImage output Fisher Information Matrix
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_fisher_grid(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *fisherpriorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation)
{
    int status = 0;
    int psf_size_x = psfImage->nx;
    int psf_size_y = psfImage->ny;
    int psf_size_semi_x = (psf_size_x-1)/2;
    int psf_size_semi_y = (psf_size_y-1)/2;
    float *fisher_matrix = (float *) fisherImage->data;
    float *fisher_matrix_prior=NULL;
    if (fisherpriorImage!=NULL)
        fisher_matrix_prior = (float *) fisherpriorImage->data;

    // 1) Project object and invert the sinogram elements
    int dim[8];
    dim[0] = 3;
    dim[1] = gridImage->dim[1];
    dim[2] = gridImage->dim[3];
    dim[3] = n_cameras;      
    dim[4] = 1;      
    dim[5] = 1;      
    dim[6] = 1;      
    dim[7] = 1;      
    nifti_image *invsinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    float *invsinogram;
    if (from_projection==0)
        {
        invsinogram = (float*) malloc(dim[1]*dim[2]*dim[3]*sizeof(float));
        invsinogramImage->data = (float *)(invsinogram);
        status = et_project(inputImage, invsinogramImage, psfImage, attenuationImage, cameras, centers, n_cameras, background, background_attenuation, 1);
        if (status)
            {
            fprintf(stderr,"'et_fisher_grid': error while calculating projection\n");
            return status;
            } 
        }
    else
        {
        invsinogramImage = inputImage;
        invsinogram = (float*)invsinogramImage->data;
        }

//    if (epsilon<=eps) epsilon=eps;
//    for (int i=0; i<invsinogramImage->nvox; i++)
//        invsinogram[i]=1/(invsinogram[i]+epsilon);

    // 2) Create (3xN) matrix of coordinates of the grid points; these will be rotated according to each camera position. 
    int n_grid_elements = fisherImage->dim[1];
    int *grid_coords  = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    for (int i=0; i<n_grid_elements*3; i++)
        grid_coords[i]=-1;
    int n=0;
    float *grid_data = (float*)gridImage->data;
    for (int x=0; x<gridImage->nx; x++) {
        for (int y=0; y<gridImage->ny; y++) {
            for (int z=0; z<gridImage->nz; z++) {
                if (grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] !=0) { //FIXME: make sure this works for non cubic images
                    n = grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] - 1;
                    if (grid_coords[n*3]!=-1)
                        return ET_ERROR_BADGRID;     // this happens if there are two elements of the grid with the same index
                    if (n > n_grid_elements) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is bigger than the numbe of elements 
                    if (n < 0) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is negative
                    grid_coords[n*3]=x;
                    grid_coords[n*3+1]=y;
                    grid_coords[n*3+2]=z;
                    }
                }
            }
        }
    // 3) For each camera position, update the FIM
    for (int i=0; i<n_grid_elements*n_grid_elements; i++)
        fisher_matrix[i]=0;
    if (fisher_matrix_prior!=NULL)
        {
        for (int i=0; i<n_grid_elements*n_grid_elements; i++)
            fisher_matrix_prior[i]=0;
        }
        
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    int *grid_coords_rotated = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    float position[3]; float position_rotated[3];
    
    for (int cam=0; cam<n_cameras; cam++)
        {
        float *invsino_data = (float *) (invsinogramImage->data) + cam*gridImage->nx*gridImage->nz ;
        // 3a) rotate the grid coordinates
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], centers[0*n_cameras+cam], centers[1*n_cameras+cam], centers[2*n_cameras+cam], XYZ_ROTATION);
        for (int n=0; n<n_grid_elements; n++)
            {
            position[0]=grid_coords[3*n]; position[1]=grid_coords[3*n+1]; position[2]=grid_coords[3*n+2]; 
            reg_mat44_mul(affineTransformation, position, position_rotated);
            grid_coords_rotated[3*n] = position_rotated[0]; grid_coords_rotated[3*n+1] = position_rotated[1]; grid_coords_rotated[3*n+2] = position_rotated[2]; //implicit rounding (NN interpolation)
            }
        // 3b) for each pair, compute the FIM element (for the current camera, then it all sums up)
        int bbox0[2];
        int bbox1[2];
        int bbox_size[2];
        int Z;
        int psf_i_x, psf_i_y, psf_j_x, psf_j_y;
        float *PSF_i; float *PSF_j;
        float fisher_element;
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords_rotated[3*i];
            int i_y = grid_coords_rotated[3*i+1];
            int i_z = grid_coords_rotated[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                int j_x = grid_coords_rotated[3*j];
                int j_y = grid_coords_rotated[3*j+1];
                int j_z = grid_coords_rotated[3*j+2];
                // 3b1) compute bounding box
                int B_x = psf_size_x - abs(j_x-i_x);     // Bounding box size
                int B_y = psf_size_y - abs(j_y-i_y);
                int p_x = max(i_x,j_x)-psf_size_semi_x;  // start of bounding box in projection space
                int p_y = max(i_y,j_y)-psf_size_semi_y;
                // 3b2) compute coordinates of the PSFs
                if (B_x>0 && B_y>0 && i_z>=0 && j_z>=0 && i_z<gridImage->nz && j_z<gridImage->nz)  // check if the PSFs overlap and if the plane if within the field of view 
                    {
                    if ((j_x >= i_x) && (j_y > i_y))
                        {
                        psf_i_x = j_x - i_x;
                        psf_i_y = j_y - i_y;
                        psf_j_x = 0;
                        psf_j_y = 0;
                        }
                    else if ((j_x < i_x) && (j_y >= i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = j_y-i_y;
                        psf_j_x = i_x-j_x;
                        psf_j_y = 0;
                        }
                    else if ((j_x <= i_x) && (j_y < i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = 0;
                        psf_j_x = i_x-j_x;
                        psf_j_y = i_y-j_y;
                        }
                    else
                        {
                        psf_i_x = j_x-i_x;
                        psf_i_y = 0;
                        psf_j_x = 0;
                        psf_j_y = i_y-j_y;
                        }         
                    // 3b3) update pointers to the PSFs (function of z)
                    PSF_i = (float*) (psfImage->data) + i_z * psf_size_x*psf_size_y; 
                    PSF_j = (float*) (psfImage->data) + j_z * psf_size_x*psf_size_y; 
                    // 3b4) update the Fisher Information matrix
                    for (int x=0; x<B_x; x++)
                        {
                        for (int y=0; y<B_y; y++)
                            {
                            // check if the point is within the projection space
                            int proj_x = p_x+x;
                            int proj_y = p_y+y;
                            if (proj_x>=0 && proj_y>=0 && proj_x<invsinogramImage->nx && proj_y<invsinogramImage->ny)
                                {
                                fisher_element = fisher_matrix[i*n_grid_elements+j] + PSF_j[(psf_j_y+y)*psf_size_x+(psf_j_x+x)]*PSF_i[(psf_i_y+y)*psf_size_x+(psf_i_x+x)]*invsino_data[proj_y*invsinogramImage->nx + proj_x];
                                //invsinogram[cam*invsinogramImage->nx*invsinogramImage->nx+proj_y*invsinogramImage->nx + proj_x];
                                fisher_matrix[i*n_grid_elements+j] = fisher_element;
                                }
                            }
                        }
                    }
                }
            }
        }

    // 3c) Fisher Information of the prior
    if (fisher_matrix_prior!=NULL)
        {
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords[3*i];
            int i_y = grid_coords[3*i+1];
            int i_z = grid_coords[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                int j_x = grid_coords[3*j];
                int j_y = grid_coords[3*j+1];
                int j_z = grid_coords[3*j+2];
                if (abs(i_x-j_x)<=1 && abs(i_y-j_y)<=1 && abs(i_z-j_z)<=1)
                    fisher_matrix_prior[i*n_grid_elements+j] = 1;
//                float dist = sqrt((i_x-j_x)^2 + (i_y-j_y)^2 + (i_z-j_z)^2);
//                fisher_matrix_prior[i*n_grid_elements+j] = dist;
                }
            }
        }

    // 4) Fill matrix (the other half)
    for (int i=0; i<n_grid_elements; i++)
        for (int j=i+1; j<n_grid_elements; j++)
            fisher_matrix[j*n_grid_elements+i] = fisher_matrix[i*n_grid_elements+j];

    if (fisher_matrix_prior!=NULL)
        for (int i=0; i<n_grid_elements; i++)
            for (int j=i+1; j<n_grid_elements; j++)
                fisher_matrix_prior[j*n_grid_elements+i] = fisher_matrix_prior[i*n_grid_elements+j];


    // Free
    free(grid_coords);
    free(affineTransformation);
    free(grid_coords_rotated);
    if (from_projection==0)
        {
        free(invsinogram);
        (void)nifti_free_extensions( invsinogramImage) ;
        free(invsinogramImage) ;
        }
    return status;
}



//! Gradient of the attenuation map. Use this in order to estimate the attenuation map from the emission data. 
/*!
  \param *gradientImage outpup gradient in voxel space. 
  \param *sinoImage input sinogram.  
  \param *activityImage input activity (estimate). 
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_gradient_attenuation(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values) 
{
    return 1;
}



//! Convolve a stack of 2D images. 
/*!
  \param *inImage input stack of images. 
  \param *outImage output convolved stack of images. 
  \param *psfImage convolution kernel. 
*/
int et_convolve(nifti_image *inImage, nifti_image *outImage, nifti_image *kernelImage)
{
    int status = 1;
    return status;
}


int et_project_partial(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *partialsumImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, nifti_image *backgroundImage, float background_attenuation, int truncate_negative_values, int do_rotate_partial)
{
return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   GPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _USE_CUDA

//! Affine transformation of nifti_image on GPU
/*!
  \param *sourceImage the source image to be transformed. 
  \param *resultImage the transformed image. 
  \param *affineTransformation the [4x4] transformed matrix. 
  \param background the background value when resampling the transformed image. 
*/
int et_affine_gpu(nifti_image *sourceImage, nifti_image *resultImage, mat44 *affineTransformation, float background)
{
    /* initialise the cuda arrays */
    cudaArray *sourceImageArray_d;
    float     *resultImageArray_d;
    float4    *positionFieldImageArray_d;
    int       *mask_d;
    
    /* Allocate arrays on the device */
    if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&resultImageArray_d, resultImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, resultImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, resultImage->dim)) return 1;

    /* Transfer data from the host to the device */
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;
    int *mask_h=(int *)malloc(resultImage->nvox*sizeof(int));
    for(int i=0; i<resultImage->nvox; i++) mask_h[i]=i;
    cudaMemcpy(mask_d, mask_h, resultImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
    free(mask_h);

    /* Apply affine */
    reg_affine_positionField_gpu(    affineTransformation,
                    resultImage,
                    &positionFieldImageArray_d);
    /* Resample the source image */
    reg_resampleSourceImage_gpu(    resultImage,
                    sourceImage,
                    &resultImageArray_d,
                    &sourceImageArray_d,
                    &positionFieldImageArray_d,
                    &mask_d,
                    resultImage->nvox,
                    background,
                    0);
    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(resultImage, &resultImageArray_d)) return 1;
    /*Free*/
    cudaCommon_free(&sourceImageArray_d);
    cudaCommon_free((void **)&resultImageArray_d);
    cudaCommon_free((void **)&mask_d);
    cudaCommon_free((void **)&positionFieldImageArray_d);

    return 0;
}


//! Rotate a nifti_image in 3D on the GPU
/*!
  \param *sourceImage the source image to be transformed. 
  \param *resultImage the rotated image. 
  \param theta_x the rotation angle around x axis in radians. 
  \param theta_y the rotation angle around y axis in radians. 
  \param theta_z the rotation angle around z axis in radians. 
  \param center_x the center of rotation along x.  
  \param center_y the center of rotation along y. 
  \param center_z the center of rotation along z. 
  \param background the background value when resampling the transformed image. 
*/
int et_rotate_gpu(nifti_image *sourceImage, nifti_image *resultImage, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, float background, int axis_order)
{
    int status = 1;
    
    //Create transformation matrix
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    et_create_rotation_matrix(affineTransformation, theta_x, theta_y, theta_z, center_x, center_y, center_z, axis_order);
    
    //Apply affine transformation
    status = et_affine_gpu(sourceImage, resultImage, affineTransformation, background);

    //Free
    free(affineTransformation);

    return status;
}


//! Projection for Emission Imaging, on GPU
/*!
  \param *activityImage the activity (or its estimate). NULL for attenuation and background activity only. 
  \param *sinoImage the photon counts in projection space. 
  \param *psfImage the depth-dependent point spread function, NULL for no PSF. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_project_gpu(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
    /* Initialise the cuda arrays */
    cudaArray *activityArray_d=NULL;               //stores input activity, makes use of fetch unit
    cudaArray *attenuationArray_d=NULL;            //stores input attenuation coefficients, makes use of fetch unit
    float     *sinoArray_d=NULL;                   //stores sinogram (output)
    float     *rotatedArray_d=NULL;                //stores activity aligned with current camera
    float     *rotatedAttenuationArray_d=NULL;     //stores attenuation coefficients aligned with current camera
    float4    *positionFieldImageArray_d=NULL;     //stores position field for rotation of activity and attenuation
    int       *mask_d=NULL;                        //binary mask that defines active voxels (typically all active)
    float     *psfArray_d=NULL;                    //stores point spread function
    float     *psfSeparatedArray_d=NULL;           //stores point spread function
    int       psf_size[3];
    int       image_size[3];
    int       separable_psf=0;

    /* Check consistency of input */
    nifti_image *referenceImage;              // this image holds information about image size and voxel size (activity or attenuation might not be defined (NULL pointers) ) 
    if (activityImage==NULL && attenuationImage==NULL){
        fprintf(stderr, "et_project_gpu: Error - define at least one between activityImage and attenuationImage. \n");
        return niftyrec_error_parameters; 
    }
    else if (attenuationImage==NULL)
        referenceImage=activityImage;
    else
        referenceImage=attenuationImage; 
    
    /* Allocate arrays on the device and transfer data to the device */
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);         
    // Activity
    if (activityImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;
        } 
        alloc_record_add(memory_record,(void*)activityArray_d,ALLOCTYPE_CUDA_ARRAY);
        if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, activityImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;
        }  
        alloc_record_add(memory_record,(void*)rotatedArray_d,ALLOCTYPE_CUDA);  
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;
        } 
    }

    // Singoram 
    if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)sinoArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, referenceImage->dim)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, referenceImage->dim)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    } 
    alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);

    // Mask 
    int *mask_h=(int *)malloc(referenceImage->nvox*sizeof(int));
    for(int i=0; i<referenceImage->nvox; i++) mask_h[i]=i;
    cudaMemcpy(mask_d, mask_h, referenceImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
    free(mask_h);
        
    /* Alloc transformation matrix */
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_HOST);

    /* Allocate and initialize kernel for DDPSF */
    if (psfImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;
        } 
        alloc_record_add(memory_record,(void*)psfArray_d,ALLOCTYPE_CUDA);
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;
        } 
        psf_size[0] = psfImage->dim[1];
        psf_size[1] = psfImage->dim[2];
        psf_size[2] = psfImage->dim[3];
        image_size[0] = referenceImage->dim[1];
        image_size[1] = referenceImage->dim[2];
        image_size[2] = referenceImage->dim[3];

        if(1)    //(psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
    }

    if (separable_psf) {
        if(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)) != cudaSuccess) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)psfSeparatedArray_d,ALLOCTYPE_CUDA);

        float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
        float psf_norm;
        for (int n=0; n<psf_size[2];n++) {
            psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
            for (int i=0;i<psf_size[0];i++) {
                psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
            }
        }
        cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice);
        free(psfSeparatedArray_h);
    }

    if (attenuationImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;
        } 
        alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
        if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)rotatedAttenuationArray_d,ALLOCTYPE_CUDA);
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_transfergpu;
        }
    }

    for(unsigned int cam=0; cam<n_cameras; cam++){
//        fprintf(stderr,"et_project_gpu: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
//        fprintf(stderr,"et_project_gpu: Center:   %f  %f  %f  \n",centers[0*n_cameras+cam], centers[1*n_cameras+cam], centers[2*n_cameras+cam]);

        // Apply affine //
        et_create_rotation_matrix(              affineTransformation, 
                                                cameras[0*n_cameras+cam], 
                                                cameras[1*n_cameras+cam], 
                                                cameras[2*n_cameras+cam], 
                                                centers[0*n_cameras+cam], 
                                                centers[1*n_cameras+cam], 
                                                centers[2*n_cameras+cam], 
                                                XYZ_ROTATION);

        reg_affine_positionField_gpu(           affineTransformation,
                                                referenceImage,
                                                &positionFieldImageArray_d);

        // Resample the activity image //
        if (activityImage != NULL)
            reg_resampleSourceImage_gpu(        activityImage,
                                                activityImage,
                                                &rotatedArray_d,
                                                &activityArray_d,
                                                &positionFieldImageArray_d,
                                                &mask_d,
                                                activityImage->nvox,
                                                background,
                                                0);

        // Resample the attenuation map //
        if (attenuationImage != NULL)
            reg_resampleSourceImage_gpu(        attenuationImage,
                                                attenuationImage,
                                                &rotatedAttenuationArray_d,
                                                &attenuationArray_d,
                                                &positionFieldImageArray_d,
                                                &mask_d,
                                                attenuationImage->nvox,
                                                background_attenuation,
                                                0);

        // Apply Depth Dependent Point Spread Function //
        if ((psfImage != NULL) && (activityImage!=NULL)) {
            if (separable_psf) {
                int status = et_convolveSeparable2D_gpu(
                                                &rotatedArray_d, 
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                if (status) {
                    alloc_record_destroy(memory_record);
                    return niftyrec_error_kernel;
                }
            }
            else
                et_convolveFFT2D_gpu(           &rotatedArray_d, 
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);
        }

        // Integrate along lines //
        if ((activityImage!=NULL) && (attenuationImage != NULL))
            et_line_integral_attenuated_gpu(    rotatedArray_d,
                                                rotatedAttenuationArray_d, 
                                                sinoArray_d,
                                                NULL,
                                                NULL,
                                                cam,
                                                referenceImage,
                                                background);
        else if ((activityImage!=NULL) && (attenuationImage == NULL))
            et_line_integral_attenuated_gpu(    rotatedArray_d,
                                                NULL,
                                                sinoArray_d,
                                                NULL,
                                                NULL,
                                                cam,
                                                referenceImage,
                                                background);
        else if ((activityImage==NULL) && (attenuationImage != NULL))
            et_line_integral_attenuated_gpu(    NULL,
                                                rotatedAttenuationArray_d, 
                                                sinoArray_d,
                                                NULL,
                                                NULL,
                                                cam,
                                                referenceImage,
                                                background);
        else {
            alloc_record_destroy(memory_record);
            return niftyrec_error_parameters;
        }
    }

    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(sinoImage, &sinoArray_d)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    }

    /* Truncate negative values: small negative values may be found due to FFT and IFFT */
    float* sino_data = (float*) sinoImage->data;
    if (truncate_negative_values) {
        for (int i=0; i<sinoImage->nvox; i++)
            if (sino_data[i] < 0)
                sino_data[i] = 0;
    }

    /*Free*/
    int status = alloc_record_destroy(memory_record); 
    return status; 
}




//! Back-projection for Emission Imaging, on GPU.
/*!
  \param *sinogramImage the data to be back-projected in projection space. 
  \param *backprojectionImage the output backprojection. 
  \param *psfImage the depth-dependent point spread function, NULL for no point spread function. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_backproject_gpu(nifti_image *sinoImage, nifti_image *backprojectionImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
    /* initialise the cuda arrays */
    cudaArray *backprojectionArray_d;
    cudaArray *attenuationArray_d;
    float     *temp_backprojection_d;
    float     *sinoArray_d;
    float     *rotatedArray_d;
    float     *rotatedAttenuationArray_d;
    float     *attenuationPlaneArray_d;
    float     *accumulatorArray_d;
    float4    *positionFieldImageArray_d;
    int       *mask_d;
    float     *psfArray_d;
    float     *psfSeparatedArray_d;
    int       psf_size[3];
    int       image_size[3];
    int       separable_psf=0;

    /* Allocate the deformation Field image */
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);  
    
    /* Allocate arrays on the device */
    cudaChannelFormatDesc backprojectionArray_d_chdesc = cudaCreateChannelDesc<float>();
    cudaExtent backprojectionArray_d_extent;
    backprojectionArray_d_extent.width  = backprojectionImage->nx;
    backprojectionArray_d_extent.height = backprojectionImage->ny;
    backprojectionArray_d_extent.depth  = backprojectionImage->nz;
    if (cudaMalloc3DArray(&backprojectionArray_d, &backprojectionArray_d_chdesc, backprojectionArray_d_extent) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)backprojectionArray_d,ALLOCTYPE_CUDA_ARRAY);

    //allocate backprojection
    cudaPitchedPtr temp_backprojection_pitched; 
    cudaExtent temp_backprojection_extent = make_cudaExtent(sizeof(float)*backprojectionImage->nx,backprojectionImage->ny,backprojectionImage->nz);     
    cudaError_t cuda_error = cudaMalloc3D(&temp_backprojection_pitched, temp_backprojection_extent);     
    if(cuda_error != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    temp_backprojection_d = (float*) temp_backprojection_pitched.ptr;
    alloc_record_add(memory_record,(void*)temp_backprojection_d,ALLOCTYPE_CUDA);

    if (attenuationImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
        if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)rotatedAttenuationArray_d,ALLOCTYPE_CUDA);
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_transfergpu; 
        }
    }

    if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)sinoArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)rotatedArray_d,ALLOCTYPE_CUDA);    
    if(cudaCommon_allocateArrayToDevice<float>(&accumulatorArray_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)accumulatorArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);

    /* Transfer data from the host to the device */
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sinoArray_d,sinoImage)) return 1;
    int *mask_h=(int *)malloc(backprojectionImage->nvox*sizeof(int)); 
    if (mask_h==NULL) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_alloccpu;
    }
    for(int i=0; i<backprojectionImage->nvox; i++) mask_h[i]=i;
    if (cudaMemcpy(mask_d, mask_h, backprojectionImage->nvox*sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_transfergpu;
    }
    free(mask_h);

    /* Alloc transformation matrix */
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    if (affineTransformation==NULL) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_alloccpu;
    }
    alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_HOST);

    /* Clear accumulator */
    et_clear_accumulator_gpu(&accumulatorArray_d,backprojectionImage );

    /* Allocate and initialize kernel for DDPSF */
    if (psfImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)psfArray_d,ALLOCTYPE_CUDA); 
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_transfergpu;
        }
        psf_size[0] = psfImage->dim[1];
        psf_size[1] = psfImage->dim[2];
        psf_size[2] = psfImage->dim[3];
        image_size[0] = backprojectionImage->dim[1];
        image_size[1] = backprojectionImage->dim[2];
        image_size[2] = backprojectionImage->dim[3];

        if(1) //(psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
    }

    if (separable_psf) {
        if (cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)psfSeparatedArray_d,ALLOCTYPE_CUDA); 
        float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
        float psf_norm;
        for (int n=0; n<psf_size[2];n++) {
            psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
            for (int i=0;i<psf_size[0];i++) {
                psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
            }
        }
        if (cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            free(psfSeparatedArray_h);
            return niftyrec_error_transfergpu;
        }
        free(psfSeparatedArray_h);
    }

    for(int cam=0; cam<n_cameras; cam++){
        // Rotate attenuation //
        if (attenuationImage != NULL) {
            et_create_rotation_matrix(          affineTransformation,
                                                cameras[0*n_cameras+cam],
                                                cameras[1*n_cameras+cam],
                                                cameras[2*n_cameras+cam],
                                                centers[0*n_cameras+cam], 
                                                centers[1*n_cameras+cam], 
                                                centers[2*n_cameras+cam], 
                                                XYZ_ROTATION);
                                                reg_affine_positionField_gpu(affineTransformation,
                                                attenuationImage,
                                                &positionFieldImageArray_d);
                                                
            reg_resampleSourceImage_gpu(        attenuationImage,
                                                attenuationImage,
                                                &rotatedAttenuationArray_d,
                                                &attenuationArray_d,
                                                &positionFieldImageArray_d,
                                                &mask_d,
                                                attenuationImage->nvox,
                                                background_attenuation,
                                                0);
        }
        // Line Backproject //
        if (attenuationImage != NULL) {
            et_line_backproject_attenuated_gpu( &sinoArray_d,
                                                &temp_backprojection_d,
                                                &rotatedAttenuationArray_d,
                                                cam,
                                                backprojectionImage);
        }
        else {
            et_line_backproject_gpu(            &sinoArray_d,
                                                &temp_backprojection_d,
                                                cam,
                                                backprojectionImage);
        }

        // Apply Depth Dependent Point Spread Function //
        if (psfImage != NULL) {
            if (separable_psf)
                et_convolveSeparable2D_gpu(     &temp_backprojection_d,
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
            else
                et_convolveFFT2D_gpu(           &temp_backprojection_d,
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
        }

        // Copy cudaArray bound memory (for rotation) //
        cudaError_t cuda_status;
        cudaExtent volumeSize = make_cudaExtent(backprojectionImage->nx, backprojectionImage->ny, backprojectionImage->nz);
        cudaMemcpy3DParms copyparms={0};
        copyparms.extent = volumeSize;
        copyparms.dstArray = backprojectionArray_d;
        copyparms.kind = cudaMemcpyDeviceToDevice; 
        copyparms.srcPtr = temp_backprojection_pitched;
        cuda_status = cudaMemcpy3D(&copyparms);
        if (cuda_status != cudaSuccess) {
            fprintf(stderr, "Error copying to texture bound memory: %s\n",cudaGetErrorString(cuda_status));
            return 1;
        }
        
        // Rotate backprojection //
        et_create_rotation_matrix(              affineTransformation,
                                                -cameras[0*n_cameras+cam],
                                                -cameras[1*n_cameras+cam],
                                                -cameras[2*n_cameras+cam],
                                                centers[0*n_cameras+cam], 
                                                centers[1*n_cameras+cam], 
                                                centers[2*n_cameras+cam], 
                                                ZYX_ROTATION);
                                                
        reg_affine_positionField_gpu(           affineTransformation,
                                                backprojectionImage,
                                                &positionFieldImageArray_d);
                                                
        reg_resampleSourceImage_gpu(            backprojectionImage,
                                                backprojectionImage,
                                                &rotatedArray_d,
                                                &backprojectionArray_d,
                                                &positionFieldImageArray_d,
                                                &mask_d,
                                                backprojectionImage->nvox,
                                                background,
                                                0);

        // Accumulate //
        et_accumulate_gpu(                      &rotatedArray_d,
                                                &accumulatorArray_d,
                                                backprojectionImage );
    }

    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(backprojectionImage, &accumulatorArray_d)) return 1; 

    /* Truncate negative values: small negative values may be found due to FFT and IFFT */
    float* accumulator_data = (float*) backprojectionImage->data;
    if (truncate_negative_values) {
        for (int i=0; i<backprojectionImage->nvox; i++)
            if (accumulator_data[i] < 0)
                accumulator_data[i] = 0;
    }

    /*Free*/
    int status = alloc_record_destroy(memory_record); 
    return status; 
}





//! Fisher Information Matrix of a grid of voxels, Emission Imaging, on GPU. 
/*!
  \param from_projection whether the input image is a projection image (1) or activity image (0)
  \param *inputImage input image: projection image or activity image. 
  \param *gridImage grid image (same size as the activity), indexes from 1 to N_points at grid locations, 0 elsewhere. 
  \param *fisherImage output Fisher Information Matrix
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_fisher_grid_gpu(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *fisherpriorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation)
{
    int status = 0;
    int psf_size_x = psfImage->nx;
    int psf_size_y = psfImage->ny;
    int psf_size_semi_x = (psf_size_x-1)/2;
    int psf_size_semi_y = (psf_size_y-1)/2;
    float *fisher_matrix = (float *) fisherImage->data;
    float *fisher_matrix_prior = NULL;
    if (fisherpriorImage!=NULL)
        fisher_matrix_prior = (float *) fisherpriorImage->data;

    // 1) Project object and invert the sinogram elements
    int dim[8];
    dim[0] = 3;
    dim[1] = gridImage->dim[1];
    dim[2] = gridImage->dim[3];
    dim[3] = n_cameras;      
    dim[4] = 1;      
    dim[5] = 1;      
    dim[6] = 1;      
    dim[7] = 1;      
    nifti_image *invsinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    float *invsinogram;
    if (from_projection==0)
        {
        invsinogram = (float*) malloc(dim[1]*dim[2]*dim[3]*sizeof(float));
        invsinogramImage->data = (float *)(invsinogram);
        status = et_project_gpu(inputImage, invsinogramImage, psfImage, attenuationImage, cameras, centers, n_cameras, background, background_attenuation, 1);
        if (status)
            {
            fprintf_verbose("'et_fisher_grid': error while calculating projection\n");
            return status;
            } 
        }
    else
        {
        invsinogramImage = inputImage;
        invsinogram = (float*)invsinogramImage->data;
        }

//    if (epsilon<=eps) epsilon=eps;
//    for (int i=0; i<invsinogramImage->nvox; i++)
//        invsinogram[i]=1/(invsinogram[i]+epsilon);

    // 2) Create (3xN) matrix of coordinates of the grid points; these will be rotated according to each camera position. 
    int n_grid_elements = fisherImage->dim[1];
    int *grid_coords  = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    for (int i=0; i<n_grid_elements*3; i++)
        grid_coords[i]=-1;
    int n=0;
    float *grid_data = (float*)gridImage->data;
    for (int x=0; x<gridImage->nx; x++) {
        for (int y=0; y<gridImage->ny; y++) {
            for (int z=0; z<gridImage->nz; z++) {
                if (grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] !=0) { //FIXME: make sure this works for non cubic images
                    n = grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] - 1;
                    if (grid_coords[n*3]!=-1)
                        return ET_ERROR_BADGRID;     // this happens if there are two elements of the grid with the same index
                    if (n > n_grid_elements) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is bigger than the numbe of elements 
                    if (n < 0) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is negative
                    grid_coords[n*3]=x;
                    grid_coords[n*3+1]=y;
                    grid_coords[n*3+2]=z;
                    }
                }
            }
        }
    // 3) For each camera position, update the FIM
    for (int i=0; i<n_grid_elements*n_grid_elements; i++)
        fisher_matrix[i]=0;
    if (fisher_matrix_prior!=NULL)
        for (int i=0; i<n_grid_elements*n_grid_elements; i++)
            fisher_matrix_prior[i]=0; 
        
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    int *grid_coords_rotated = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    float position[3]; float position_rotated[3];
    
    for (int cam=0; cam<n_cameras; cam++)
        {
        float *invsino_data = (float *) (invsinogramImage->data) + cam*gridImage->nx*gridImage->nz ;
        // 3a) rotate the grid coordinates
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], centers[0*n_cameras+cam], centers[1*n_cameras+cam], centers[2*n_cameras+cam], XYZ_ROTATION);
        for (int n=0; n<n_grid_elements; n++)
            {
            position[0]=grid_coords[3*n]; position[1]=grid_coords[3*n+1]; position[2]=grid_coords[3*n+2]; 
            reg_mat44_mul(affineTransformation, position, position_rotated);
            grid_coords_rotated[3*n] = position_rotated[0]; grid_coords_rotated[3*n+1] = position_rotated[1]; grid_coords_rotated[3*n+2] = position_rotated[2]; //implicit rounding (NN interpolation)
            }
        // 3b) for each pair, compute the FIM element (for the current camera, then it all sums up)
        int bbox0[2];
        int bbox1[2];
        int bbox_size[2];
        int Z;
        int psf_i_x, psf_i_y, psf_j_x, psf_j_y;
        float *PSF_i; float *PSF_j;
        float fisher_element;
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords_rotated[3*i];
            int i_y = grid_coords_rotated[3*i+1];
            int i_z = grid_coords_rotated[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                int j_x = grid_coords_rotated[3*j];
                int j_y = grid_coords_rotated[3*j+1];
                int j_z = grid_coords_rotated[3*j+2];
                // 3b1) compute bounding box
                int B_x = psf_size_x - abs(j_x-i_x);     // Bounding box size
                int B_y = psf_size_y - abs(j_y-i_y);
                int p_x = max(i_x,j_x)-psf_size_semi_x;  // start of bounding box in projection space
                int p_y = max(i_y,j_y)-psf_size_semi_y;
                // 3b2) compute coordinates of the PSFs
                if (B_x>0 && B_y>0 && i_z>=0 && j_z>=0 && i_z<gridImage->nz && j_z<gridImage->nz)  // check if the PSFs overlap and if the plane if within the field of view 
                    {
                    if ((j_x >= i_x) && (j_y > i_y))
                        {
                        psf_i_x = j_x - i_x;
                        psf_i_y = j_y - i_y;
                        psf_j_x = 0;
                        psf_j_y = 0;
                        }
                    else if ((j_x < i_x) && (j_y >= i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = j_y-i_y;
                        psf_j_x = i_x-j_x;
                        psf_j_y = 0;
                        }
                    else if ((j_x <= i_x) && (j_y < i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = 0;
                        psf_j_x = i_x-j_x;
                        psf_j_y = i_y-j_y;
                        }
                    else
                        {
                        psf_i_x = j_x-i_x;
                        psf_i_y = 0;
                        psf_j_x = 0;
                        psf_j_y = i_y-j_y;
                        }         
                    // 3b3) update pointers to the PSFs (function of z)
                    PSF_i = (float*) (psfImage->data) + i_z * psf_size_x*psf_size_y; 
                    PSF_j = (float*) (psfImage->data) + j_z * psf_size_x*psf_size_y; 
                    // 3b4) update the Fisher Information matrix
                    for (int x=0; x<B_x; x++)
                        {
                        for (int y=0; y<B_y; y++)
                            {
                            // check if the point is within the projection space
                            int proj_x = p_x+x;
                            int proj_y = p_y+y;
                            if (proj_x>=0 && proj_y>=0 && proj_x<invsinogramImage->nx && proj_y<invsinogramImage->ny)
                                {
                                fisher_element = fisher_matrix[i*n_grid_elements+j] + PSF_j[(psf_j_y+y)*psf_size_x+(psf_j_x+x)]*PSF_i[(psf_i_y+y)*psf_size_x+(psf_i_x+x)]*invsino_data[proj_y*invsinogramImage->nx + proj_x];
                                //invsinogram[cam*invsinogramImage->nx*invsinogramImage->nx+proj_y*invsinogramImage->nx + proj_x];
                                fisher_matrix[i*n_grid_elements+j] = fisher_element;
                                }
                            }
                        }
                    }
                }
            }
        }

    // 3c) Fisher Information of the prior
    if (fisher_matrix_prior!=NULL)
        {
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords[3*i];
            int i_y = grid_coords[3*i+1];
            int i_z = grid_coords[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                if (i!=j)
                    {
                    int j_x = grid_coords[3*j];
                    int j_y = grid_coords[3*j+1];
                    int j_z = grid_coords[3*j+2];
                    if (abs(i_x-j_x)<=3 && abs(i_y-j_y)<=3 && abs(i_z-j_z)<=3)
                        fisher_matrix_prior[i*n_grid_elements+j] = 1;
                    }
//                float dist = sqrt((i_x-j_x)^2 + (i_y-j_y)^2 + (i_z-j_z)^2);
//                fisher_matrix_prior[i*n_grid_elements+j] = dist;
                }
            }
        }

    // 4) Fill matrix (the other half)
    for (int i=0; i<n_grid_elements; i++)
        for (int j=i+1; j<n_grid_elements; j++)
            fisher_matrix[j*n_grid_elements+i] = fisher_matrix[i*n_grid_elements+j];

    if (fisher_matrix_prior!=NULL)
        for (int i=0; i<n_grid_elements; i++)
            for (int j=i+1; j<n_grid_elements; j++)
                fisher_matrix_prior[j*n_grid_elements+i] = fisher_matrix_prior[i*n_grid_elements+j];

    // Free
    free(grid_coords);
    free(affineTransformation);
    free(grid_coords_rotated);
    if (from_projection==0)
        {
        free(invsinogram);
        (void)nifti_free_extensions( invsinogramImage) ;
        free(invsinogramImage) ;
        }
    return status;
} 



//! Gradient of the attenuation map. Use this in order to estimate the attenuation map from the emission data. GPU version. 
/*!
  \param *gradientImage outpup gradient in voxel space. 
  \param *sinoImage input sinogram.  
  \param *activityImage input activity (estimate). 
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_gradient_attenuation_gpu(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, float background_attenuation, int truncate_negative_values) 
{
    /* initialise the cuda arrays */
    cudaArray *activityArray_d;               //stores input activity, makes use of fetch unit
    cudaArray *backprojectionArray_d;
    cudaArray *attenuationArray_d;
    float     *temp_backprojection_d;
    float     *sinoArray_d;
    float     *rotatedArray_d;
    float     *rotatedAttenuationArray_d;
        float     *attenuationPlaneArray_d;
    float     *gradientArray_d;
    float4    *positionFieldImageArray_d;
    int       *mask_d;
        float     *psfArray_d;
        float     *psfSeparatedArray_d;
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

    /* Allocate the deformation Field image */
    nifti_image *positionFieldImage = nifti_copy_nim_info(gradientImage);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[1]=positionFieldImage->nx = gradientImage->nx;
    positionFieldImage->dim[2]=positionFieldImage->ny = gradientImage->ny;
    positionFieldImage->dim[3]=positionFieldImage->nz = gradientImage->nz;
    positionFieldImage->dim[4]=positionFieldImage->nt = 1; positionFieldImage->pixdim[4]=positionFieldImage->dt = 1.0;
    positionFieldImage->dim[5]=positionFieldImage->nu = 3; positionFieldImage->pixdim[5]=positionFieldImage->du = 1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv = 1; positionFieldImage->pixdim[6]=positionFieldImage->dv = 1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw = 1; positionFieldImage->pixdim[7]=positionFieldImage->dw = 1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    positionFieldImage->nbyper = sizeof(float);
    positionFieldImage->data=NULL;
    
    /* Allocate arrays on the device */
        cudaChannelFormatDesc backprojectionArray_d_chdesc = cudaCreateChannelDesc<float>();
        cudaExtent backprojectionArray_d_extent;
        backprojectionArray_d_extent.width  = gradientImage->nx;
        backprojectionArray_d_extent.height = gradientImage->ny;
        backprojectionArray_d_extent.depth  = gradientImage->nz;
        cudaError_t cuda_status1 = cudaMalloc3DArray(&backprojectionArray_d, &backprojectionArray_d_chdesc, backprojectionArray_d_extent);

    if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) return 1;

        cudaPitchedPtr temp_backprojection_pitched; 
        cudaExtent temp_backprojection_extent = make_cudaExtent(sizeof(float)*gradientImage->nx,gradientImage->ny,gradientImage->nz);     
        cudaError_t cuda_error = cudaMalloc3D(&temp_backprojection_pitched, temp_backprojection_extent);     
        if(cuda_error != cudaSuccess) {
            return niftyrec_error_allocgpu;}
        temp_backprojection_d = (float*) temp_backprojection_pitched.ptr;

        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) return 1;
        if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) return 1;    
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) return 1;

    if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, gradientImage->dim)) return 1;    
    if(cudaCommon_allocateArrayToDevice<float>(&gradientArray_d, gradientImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, gradientImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, gradientImage->dim)) return 1;

    /* Transfer data from the host to the device */
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) return 1;
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sinoArray_d,sinoImage)) return 1;
    int *mask_h=(int *)malloc(gradientImage->nvox*sizeof(int));
    for(int i=0; i<gradientImage->nvox; i++) mask_h[i]=i;
    cudaMemcpy(mask_d, mask_h, gradientImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
    free(mask_h);

    /* Alloc transformation matrix */
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

    /* Clear gradient */
    et_clear_accumulator_gpu(        &gradientArray_d,
                        gradientImage );
        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = gradientImage->dim[1];
            image_size[1] = gradientImage->dim[2];
            image_size[2] = gradientImage->dim[3];

            if(1) //(psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice);
            free(psfSeparatedArray_h);
            }

    for(int cam=0; cam<n_cameras; cam++){
                // Rotate attenuation and activity //
                et_create_rotation_matrix(      affineTransformation,
                        cameras[0*n_cameras+cam],
                        cameras[1*n_cameras+cam],
                        cameras[2*n_cameras+cam],
                        centers[0*n_cameras+cam],
                        centers[1*n_cameras+cam],
                        centers[2*n_cameras+cam],
                        XYZ_ROTATION);
                reg_affine_positionField_gpu(   affineTransformation,
                        attenuationImage,
                        &positionFieldImageArray_d);
                reg_resampleSourceImage_gpu(    attenuationImage,
                        attenuationImage,
                        &rotatedAttenuationArray_d,
                        &attenuationArray_d,
                        &positionFieldImageArray_d,
                        &mask_d,
                        attenuationImage->nvox,
                        background_attenuation,
                        0);
                reg_resampleSourceImage_gpu(    activityImage,
                        activityImage,
                        &rotatedArray_d,
                        &activityArray_d,
                        &positionFieldImageArray_d,
                        &mask_d,
                        activityImage->nvox,
                        background,
                        0);

        // Line Backproject, compute gradient of the likelihood with respect of the attenuation coefficients //
                et_attenuation_gradient_gpu(    &rotatedArray_d, 
                                                &sinoArray_d, 
                        &temp_backprojection_d, 
                        &rotatedAttenuationArray_d, 
                                                cam, 
                        gradientImage); 


                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D_gpu( &temp_backprojection_d,
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
                    else
                        et_convolveFFT2D_gpu(   &temp_backprojection_d,
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
                    }

                // Copy to texture bound memory (for rotation) //
                cudaError_t cuda_status;
                cudaExtent volumeSize = make_cudaExtent(gradientImage->nx, gradientImage->ny, gradientImage->nz);
                cudaMemcpy3DParms copyparms={0};

                copyparms.extent = volumeSize;
                copyparms.dstArray = backprojectionArray_d;
                copyparms.kind = cudaMemcpyDeviceToDevice; 
                copyparms.srcPtr = temp_backprojection_pitched;
                cuda_status = cudaMemcpy3D(&copyparms);
                if (cuda_status != cudaSuccess)
                        {
                        fprintf(stderr, "Error copying to texture bound memory: %s\n",cudaGetErrorString(cuda_status));
                        return 1;
                        }

        // Rotate backprojection //
        et_create_rotation_matrix(    affineTransformation,
                        -cameras[0*n_cameras+cam],
                        -cameras[1*n_cameras+cam],
                        -cameras[2*n_cameras+cam],
                        centers[0*n_cameras+cam],
                        centers[1*n_cameras+cam],
                        centers[2*n_cameras+cam],
                        ZYX_ROTATION);
        reg_affine_positionField_gpu(    affineTransformation,
                        gradientImage,
                        &positionFieldImageArray_d);
        reg_resampleSourceImage_gpu(    gradientImage,
                        gradientImage,
                        &rotatedArray_d,
                        &backprojectionArray_d,
                        &positionFieldImageArray_d,
                        &mask_d,
                        gradientImage->nvox,
                        background,
                        0);

        // Accumulate //
        et_accumulate_gpu(        &rotatedArray_d,
                        &gradientArray_d,
                        gradientImage );
    }

    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(gradientImage, &gradientArray_d)) return 1; 

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* gradient_data = (float*) gradientImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<gradientImage->nvox; i++)
                if (gradient_data[i] < 0)
                    gradient_data[i] = 0;
            }

    /*Free*/
    cudaCommon_free(&activityArray_d);
        cudaCommon_free(&backprojectionArray_d);

        cudaCommon_free(&attenuationArray_d);
        cudaCommon_free((void **)&rotatedAttenuationArray_d);

    cudaCommon_free((void **)&rotatedArray_d);
    cudaCommon_free((void **)&sinoArray_d);
    cudaCommon_free((void **)&gradientArray_d);
    cudaCommon_free((void **)&mask_d);
    cudaCommon_free((void **)&positionFieldImageArray_d);
    cudaCommon_free((void **)&temp_backprojection_d);
        if (psfImage != NULL)
            {
            cudaCommon_free((void **)&psfArray_d);
            if (separable_psf)
                cudaCommon_free((void **)&psfSeparatedArray_d);
            }
    nifti_image_free(positionFieldImage);
    free(affineTransformation);

    return 0;
}



//! Convolve a stack of 2D images on the GPU. 
/*!
  \param *inImage input stack of images. 
  \param *outImage output convolved stack of images. 
  \param *psfImage convolution kernel. 
*/
int et_convolve_gpu(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage)
{
    int status = 1;

    int data_size[3];
    int kernel_size[3];

    float *d_Input;
    float *d_Kernel;
    float *d_Result;

    fprintf_verbose("et_convolve_gpu - Image size: %d %d %d\n",inImage->dim[1],inImage->dim[2],inImage->dim[3]);

    fprintf_verbose("et_convolve_gpu - Allocating memory...\n");
    if(cudaCommon_allocateArrayToDevice<float>(&d_Input, inImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&d_Kernel, psfImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&d_Result, outImage->dim)) return 1;

    fprintf_verbose("et_convolve_gpu - Uploading to GPU device...\n");
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&d_Input, inImage)) return 1;
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&d_Kernel, psfImage)) return 1;
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&d_Result, outImage)) return 1;

    fprintf_verbose("et_convolve_gpu - Performing 2D convolution...\n");
    data_size[0] = inImage->dim[0];
    data_size[1] = inImage->dim[1];
    data_size[2] = inImage->dim[2];
    kernel_size[0] = psfImage->dim[0];
    kernel_size[1] = psfImage->dim[1];
    kernel_size[2] = psfImage->dim[2];
    
    et_convolveFFT2D_gpu(&d_Input, data_size, &d_Kernel, kernel_size, &d_Result);

    fprintf_verbose("et_convolve_gpu - Reading back GPU FFT results...\n");
    //CUDA_SAFE_CALL( cudaMemcpy(h_Result, d_PaddedData, fftH * fftW * sizeof(float), cudaMemcpyDeviceToHost) );
    if(cudaCommon_transferFromDeviceToNifti(outImage, &d_Result)) return 1; 

    fprintf_verbose("et_convolve_gpu - Freeign GPU memory...\n");    
    CUDA_SAFE_CALL( cudaFree(d_Input) );
    CUDA_SAFE_CALL( cudaFree(d_Kernel) );
    CUDA_SAFE_CALL( cudaFree(d_Result) );
    fprintf_verbose("et_convolve_gpu - Freeign GPU memory...\n");   

    status = 0;
    return status;
}



//! List NVIDIA CUDA compatible GPU's installed in the system. 
/*!
  \param *device_count_out output, number of installed GPUs. 
  \param *devices outoput, GPU devices compute capability and ID's. See et_list_gpus_mex.  
*/
int et_list_gpus(int *device_count_out, int *devices)
{
    /* Initialise the cuda card */
    int status = 1;
    struct cudaDeviceProp deviceProp;
    int device_count = 0;
    int multiprocessors = 0;
    int clock = 0;
    int gflops = 0;
    int globalmem = 0;

    cudaGetDeviceCount( &device_count );

    int device_count_max = device_count;
    if (device_count_max > MAX_DEVICES)
        device_count_max = MAX_DEVICES;

    for (int dev=0; dev<device_count_max; dev++)
        {
        cudaGetDeviceProperties(&deviceProp, dev);
        multiprocessors = deviceProp.multiProcessorCount;
        clock = deviceProp.clockRate;
        gflops = multiprocessors * clock;
        globalmem = (int)floor(deviceProp.totalGlobalMem/1000000.0);
        //fprintf_verbose("\nDevice %d: %d MP, %d GHz, %d GFlops, %d Mb",dev,multiprocessors,clock,gflops,globalmem);
        devices[SIZE_OF_INFO*dev+0] = dev;
        devices[SIZE_OF_INFO*dev+1] = gflops;
        devices[SIZE_OF_INFO*dev+2] = multiprocessors;
        devices[SIZE_OF_INFO*dev+3] = clock;
        devices[SIZE_OF_INFO*dev+4] = globalmem;
        }
    device_count_out[0] = device_count;
    status = 0;
    return status;
}


//! Set GPU to be used by NiftyRec. 
/*!
  \param id GPU ID. 
*/
int et_set_gpu(int id)
{
    int status = 1;
    struct cudaDeviceProp deviceProp;

    cudaSetDevice( id );
    cudaGetDeviceProperties(&deviceProp, id );
    if (deviceProp.major < 1)
        {
        printf("ERROR - The specified graphical card does not exist.\n");
        status = 1;
    }
    else
        status = 0;
    return status;
}






/*
int et_joint_histogram_gpu(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B)
{
    float *matrix_A_d;
    float *matrix_B_d;
    int *joint_histogram_d;
    
    // Allocate arrays on device 
    cudaMalloc((void **)&matrix_A_d,  matrix_A_Image->nvox*sizeof(float));
    cudaMalloc((void **)&matrix_B_d,  matrix_B_Image->nvox*sizeof(float));
    cudaMalloc((void **)&joint_histogram_d,  joint_histogram_Image->nvox*sizeof(int));
        cudaMemset((void*)joint_histogram_d,0,joint_histogram_Image->nvox*sizeof(int));
        
    // Transfer data from the host to the device 
    cudaMemcpy(matrix_A_d, matrix_A_Image->data, matrix_A_Image->nvox*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(matrix_B_d, matrix_B_Image->data, matrix_B_Image->nvox*sizeof(float), cudaMemcpyHostToDevice);

    // Compute joint histogram 
    et_joint_histogram_gpu(        &matrix_A_d,
                    &matrix_B_d,
                    &joint_histogram_d,
                    matrix_A_Image->nvox,
                    joint_histogram_Image->nx,
                    min_A,
                    max_A,
                    min_B,
                    max_B );
    
    // Transfer joint histogram back to host 
    fprintf_verbose( "\nJH: %d %d %d %f %f %f %f",joint_histogram_Image->nx, joint_histogram_Image->nvox, matrix_A_Image->nvox, min_A, max_A, min_B, max_B);
    cudaMemcpy(joint_histogram_Image->data, joint_histogram_d, joint_histogram_Image->nvox*sizeof(int), cudaMemcpyDeviceToHost);
    
    // Free arrays on device 
    cudaCommon_free((void **)&matrix_A_d);
    cudaCommon_free((void **)&matrix_B_d);
    cudaCommon_free((void **)&joint_histogram_d);
    
    return 0;
}
*/



//! Projection for Emission Imaging, on GPU
/*!
  \param *activityImage the activity (or its estimate). NULL for attenuation and background activity only. 
  \param *sinoImage the photon counts in projection space. 
  \param *psfImage the depth-dependent point spread function, NULL for no PSF. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param *centers [n_camerasx3] array of center of the cameras in voxels. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_project_partial_gpu(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *partialsumImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, float *centers, int n_cameras, float background, nifti_image *backgroundImage, float background_attenuation, int truncate_negative_values, int do_rotate_partial)
{
    /* initialise the cuda arrays */
    cudaArray *activityArray_d=NULL;               //stores input activity, makes use of fetch unit
        cudaArray *attenuationArray_d=NULL;            //stores input attenuation coefficients, makes use of fetch unit
    float     *sinoArray_d=NULL;                   //stores sinogram (output)
    float     *rotatedArray_d=NULL;                //stores activity aligned with current camera
        float     *rotatedAttenuationArray_d=NULL;     //stores attenuation coefficients aligned with current camera
    float4    *positionFieldImageArray_d=NULL;     //stores position field for rotation of activity and attenuation
    int       *mask_d=NULL;                        //binary mask that defines active voxels (typically all active)
        float     *psfArray_d=NULL;                    //stores point spread function
        float     *psfSeparatedArray_d=NULL;           //stores point spread function
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

        /* Check consistency of input */
        nifti_image *referenceImage;              // this image holds information about image size and voxel size (activity or attenuation might not be defined (NULL pointers) ) 
        if (activityImage==NULL && attenuationImage==NULL)
            {
            fprintf(stderr, "et_project_partial_gpu: Error - define at least one between activityImage and attenuationImage. \n");
            return niftyrec_error_parameters; 
            }
        else if (attenuationImage==NULL)
            referenceImage=activityImage;
        else
            referenceImage=attenuationImage; 
    
    /* Allocate arrays on the device and transfer data to the device */

        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);         
        // Activity
        if (activityImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;} 
            alloc_record_add(memory_record,(void*)activityArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, activityImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}  
            alloc_record_add(memory_record,(void*)rotatedArray_d,ALLOCTYPE_CUDA);  
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_transfergpu;} 
            }
        // Partial sum
        //partialsumArray_d: stores the partial integrals. Pitched as it is then transfered to CUDA Array for fast rotation
    float     *partialsumArray_d=NULL;             
        cudaPitchedPtr partialsumArray_pitched; 
        cudaExtent partialsumArray_extent = make_cudaExtent(sizeof(float)*referenceImage->nx,referenceImage->ny,referenceImage->nz);     
        cudaError_t cuda_error = cudaMalloc3D(&partialsumArray_pitched, partialsumArray_extent);     
        if(cuda_error != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        partialsumArray_d = (float*) partialsumArray_pitched.ptr;
        alloc_record_add(memory_record,(void*)partialsumArray_d,ALLOCTYPE_CUDA);

        //temp_partialsumArray_d: CUDA array, temporary for fast rotation 
    cudaArray *temp_partialsumArray_d=NULL;
        cudaChannelFormatDesc temp_partialsumArray_d_chdesc = cudaCreateChannelDesc<float>();
        cudaExtent temp_partialsumArray_d_extent;
        temp_partialsumArray_d_extent.width  = referenceImage->nx;
        temp_partialsumArray_d_extent.height = referenceImage->ny;
        temp_partialsumArray_d_extent.depth  = referenceImage->nz;
        if (cudaMalloc3DArray(&temp_partialsumArray_d, &temp_partialsumArray_d_chdesc, temp_partialsumArray_d_extent) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)temp_partialsumArray_d,ALLOCTYPE_CUDA_ARRAY);

        //partialsumRotatedArray_d: partial sums rotated to the input image frame
    float     *partialsumRotatedArray_d;  
        if(cudaCommon_allocateArrayToDevice<float>(&partialsumRotatedArray_d, referenceImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)partialsumRotatedArray_d,ALLOCTYPE_CUDA);

        // Singoram 
    if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)sinoArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;} 
        alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);

    // Mask 
    int *mask_h=(int *)malloc(referenceImage->nvox*sizeof(int));
    for(int i=0; i<referenceImage->nvox; i++) mask_h[i]=i;
    cudaMemcpy(mask_d, mask_h, referenceImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
    free(mask_h);
        
    /* Alloc transformation matrix */
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_HOST);

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;} 
            alloc_record_add(memory_record,(void*)psfArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_transfergpu;} 
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = referenceImage->dim[1];
            image_size[1] = referenceImage->dim[2];
            image_size[2] = referenceImage->dim[3];

            if(1) //(psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            if(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)) != cudaSuccess) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)psfSeparatedArray_d,ALLOCTYPE_CUDA);

            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice);
            free(psfSeparatedArray_h);
            }

        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) {
                 alloc_record_destroy(memory_record); 
                 return niftyrec_error_transfergpu;} 
            alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) {
                 alloc_record_destroy(memory_record); 
                 return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)rotatedAttenuationArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) {
                 alloc_record_destroy(memory_record);
                 return niftyrec_error_transfergpu;}
            }

        float *backgroundArray_d=NULL; 
        if (backgroundImage != NULL)
            {
            if(cudaMalloc((void **)&backgroundArray_d, backgroundImage->nx*backgroundImage->ny*sizeof(float)) != cudaSuccess) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)backgroundArray_d,ALLOCTYPE_CUDA);

//            if(cudaCommon_allocateArrayToDevice<float>(&backgroundArray_d, backgroundImage->dim)) {
//                 alloc_record_destroy(memory_record); 
//                 return niftyrec_error_transfergpu;} 
//            alloc_record_add(memory_record,(void*)backgroundArray_d,ALLOCTYPE_CUDA);

            cudaMemcpy((float*)backgroundArray_d,(float*)backgroundImage->data, backgroundImage->nx*backgroundImage->ny*sizeof(float), cudaMemcpyHostToDevice); 
//            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&backgroundArray_d,backgroundImage)) {
//                 alloc_record_destroy(memory_record);
//                 return niftyrec_error_transfergpu;}
            }


    for(unsigned int cam=0; cam<n_cameras; cam++){
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
        // Apply affine //
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], centers[0*n_cameras+cam], centers[1*n_cameras+cam], centers[2*n_cameras+cam], XYZ_ROTATION);
        reg_affine_positionField_gpu(    affineTransformation,
                        referenceImage,
                        &positionFieldImageArray_d);

        // Resample the activity image //
                if (activityImage != NULL)
            reg_resampleSourceImage_gpu(activityImage,
                        activityImage,
                        &rotatedArray_d,
                        &activityArray_d,
                        &positionFieldImageArray_d,
                        &mask_d,
                        activityImage->nvox,
                        background,
                        0);

                // Resample the attenuation map //
                if (attenuationImage != NULL)
            reg_resampleSourceImage_gpu(attenuationImage,
                        attenuationImage,
                        &rotatedAttenuationArray_d,
                        &attenuationArray_d,
                        &positionFieldImageArray_d,
                        &mask_d,
                        attenuationImage->nvox,
                        background_attenuation,
                        0);

                // Apply Depth Dependent Point Spread Function //
                if ((psfImage != NULL) && (activityImage!=NULL))
                    {
                    if (separable_psf)
                        {
                        int status = et_convolveSeparable2D_gpu(
                                                &rotatedArray_d, 
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                        if (status)
                            {
                            alloc_record_destroy(memory_record);
                            return niftyrec_error_kernel;
                            }
                        }
                    else
                        et_convolveFFT2D_gpu(   &rotatedArray_d, 
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                    }

        // Integrate along lines //
                if ((activityImage!=NULL) && (attenuationImage != NULL))
                    et_line_integral_attenuated_gpu(    rotatedArray_d,
                        rotatedAttenuationArray_d, 
                        sinoArray_d,
                                                backgroundArray_d, 
                                                partialsumArray_d,
                        cam,
                        referenceImage,
                                                background);
                else if ((activityImage!=NULL) && (attenuationImage == NULL))
            et_line_integral_attenuated_gpu(    rotatedArray_d,
                                                NULL,
                        sinoArray_d,
                                                backgroundArray_d, 
                                                partialsumArray_d,
                        cam,
                        referenceImage,
                                                background);
                else if ((activityImage==NULL) && (attenuationImage != NULL))
                    et_line_integral_attenuated_gpu(NULL,
                        rotatedAttenuationArray_d, 
                        sinoArray_d,
                                                backgroundArray_d, 
                                                partialsumArray_d,
                        cam,
                        referenceImage,
                                                background);
                else 
                    {
                    alloc_record_destroy(memory_record);
                    return niftyrec_error_parameters;
                    }
                if (do_rotate_partial)
                {
                //Rotate partial integrals back to image frame
                cudaError_t cuda_status;
                cudaExtent volumeSize = make_cudaExtent(referenceImage->nx, referenceImage->ny, referenceImage->nz);
                cudaMemcpy3DParms copyparms={0};
                copyparms.extent = volumeSize;
                copyparms.dstArray = temp_partialsumArray_d; 
                copyparms.kind = cudaMemcpyDeviceToDevice; 
                copyparms.srcPtr = partialsumArray_pitched;
                cuda_status = cudaMemcpy3D(&copyparms);
                if (cuda_status != cudaSuccess)
                        {
                        fprintf(stderr, "Error copying to texture bound memory: %s\n",cudaGetErrorString(cuda_status));
                        return 1;
                        }


        et_create_rotation_matrix(    affineTransformation,
                                                -cameras[0*n_cameras+cam],
                                                -cameras[1*n_cameras+cam],
                                                -cameras[2*n_cameras+cam],
                                                centers[0*n_cameras+cam],
                                                centers[1*n_cameras+cam],
                                                centers[2*n_cameras+cam],
                                                ZYX_ROTATION);
        reg_affine_positionField_gpu(    affineTransformation,
                                                referenceImage,
                                                &positionFieldImageArray_d);
        reg_resampleSourceImage_gpu(    referenceImage,
                                                referenceImage,
                                                &partialsumRotatedArray_d,
                                                &temp_partialsumArray_d,
                                                &positionFieldImageArray_d,
                                                &mask_d,
                                                referenceImage->nvox,
                                                background,
                                                0);

                if (cudaMemcpy(((float*) partialsumImage->data)+referenceImage->nx*referenceImage->ny*referenceImage->nz*cam, partialsumRotatedArray_d, referenceImage->nx*referenceImage->ny*referenceImage->nz*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
                    alloc_record_destroy(memory_record); 
                    return niftyrec_error_transfergpu;}
                }
                else
                {
                if (cudaMemcpy(((float*) partialsumImage->data)+referenceImage->nx*referenceImage->ny*referenceImage->nz*cam, partialsumArray_d, referenceImage->nx*referenceImage->ny*referenceImage->nz*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
                    alloc_record_destroy(memory_record); 
                    return niftyrec_error_transfergpu;}
                }
    }

    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(sinoImage, &sinoArray_d)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<sinoImage->nvox; i++)
                if (sino_data[i] < 0)
                    sino_data[i] = 0;
            }

    /*Free*/
        return alloc_record_destroy(memory_record); 
}



//! Reset GPU
int et_reset_gpu()
{
    #ifdef _WIN32
    return niftyrec_error_nogpubuilt;  
    #else
    cudaDeviceReset();
    return niftyrec_success;
    #endif
}

#endif




#ifdef _USE_CUDA
//***************************************************************************************
/* PET */
unsigned int PET_project_compressed_gpu(nifti_image *activityImage, nifti_image *attenuationImage, float *projection, 
               int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, 
               unsigned int N_axial, unsigned int N_azimuthal, 
               float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
               unsigned int N_samples, float sample_step, float background, float background_attenuation, 
               unsigned int truncate_negative_values, unsigned int direction, unsigned int block_size, float *time_profiling)
{
    // performance parameters: default values ('direction' only affects performance, identical results with all 6 values; 'block_size' only affects performance - it's the number of parallel threads per block)
    if (block_size==0)
        block_size = 512; 

    // N_locations:                                           number of locations of active interaction (compressed projection data structure)
    // N_axial:                                               number of detector planes, axial 
    // N_azimuthal:                                           number of detector planes, azimuthal  

    // N_samples:                                             number of samples in direction perpendicular to the detector plane
    // sample_step:                                           distance between sample points in direction perpendicular to detector plane 
    // background:                                            background activity (in each voxel of the background - used when resampling) 
    // background_attenuation:                                background attenuation map (in each voxel of the background - used when resampling)
    // truncate_negative_values:                              flag: if 1 then set negative values in the resulting projection to 0 

    struct timeval t1, t2; 
    double T0_transfer_to_gpu, T1_alloc, T2_rotate, T3_resample, T4_integral, T5_transfer_to_host, T6_free;
    T0_transfer_to_gpu = 0.0; 
    T1_alloc = 0.0; 
    T2_rotate = 0.0; 
    T3_resample = 0.0;
    T4_integral = 0.0; 
    T5_transfer_to_host = 0.0;
    T6_free = 0.0; 
    

    /* define variables - note that *_d indicates variables that reside on the cuda (d)evice, they are of two types: cudaArray (for texture fetching) or not cudaArray */
    cudaArray      *activityArray_d             = NULL;     //stores input activity, makes use of GPU texture fetch unit
    cudaArray      *attenuationArray_d          = NULL;     //stores input attenuation coefficients, makes use of GPU texture fetch unit
    float          *resampledActivityArray_d    = NULL;     //temporarily stores activity aligned with one of the detector plane  
    float          *resampledAttenuationArray_d = NULL;     //temporarily stores attenuation coefficients aligned with one of the detector planes
    float4         *positionFieldImageArray_d   = NULL;     //temporarily stores position field for rotation of activity and attenuation 
    int            *mask_d                      = NULL;     //binary mask that defines active voxels (typically all active)
    unsigned short *locationsArray_d            = NULL;     //stores the locations array (compressed projection data structure)
    float          *projectionArray_d           = NULL;     //stores the projection data (output)

    /* Allocate arrays on the device and transfer data to the device */
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);  //memory handling: keep a record of the allocated memory 

    // -- Activity --
    gettimeofday(&t1, NULL);
    if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough memory on GPU (%lu bytes) for 'activityArray_d'. \n",activityImage->ndim*sizeof(float));
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)activityArray_d,ALLOCTYPE_CUDA_ARRAY);
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 
    
    
    gettimeofday(&t1, NULL); 
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Unable to copy %lu bytes of memory to GPU: from 'activityImage->data' to 'activityArray_d'. \n",activityImage->ndim*sizeof(float));
        return niftyrec_error_transfergpu;
    }
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    
    // -- Resampled activity -- 
    //  1) define sform matrix: 
    gettimeofday(&t1, NULL);
    mat44 *scale_m       = (mat44 *)calloc(1,sizeof(mat44));
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44));
    mat44 *m             = (mat44 *)calloc(1,sizeof(mat44));
    alloc_record_add(memory_record,(void*)scale_m,ALLOCTYPE_HOST);
    alloc_record_add(memory_record,(void*)translation_m,ALLOCTYPE_HOST);
    alloc_record_add(memory_record,(void*)m,ALLOCTYPE_HOST);
    
    // 2) define resample activity nifti image
    int dim[8]; 
    dim[0]    = 3; 
    if (direction==1) {
        dim[1] = N_u;                dim[2] = N_v;              dim[3] = N_samples; 
        et_create_scale_matrix(         scale_m,         size_u/N_u,         size_v/N_v,        sample_step); 
        et_create_translation_matrix(   translation_m,   -0.5*size_u,        -0.5*size_v,       -0.5 * sample_step * N_samples);
    }
    else if (direction==2) {
        dim[1] = N_v;                dim[2] = N_u;              dim[3] = N_samples; 
        et_create_scale_matrix(         scale_m,         size_v/N_v,         size_u/N_u,        sample_step); 
        et_create_translation_matrix(   translation_m,   -0.5*size_v,        -0.5*size_u,       -0.5 * sample_step * N_samples);
    }
    else if (direction==3) {
        dim[1] = N_samples;          dim[2] = N_u;              dim[3] = N_v; 
        et_create_scale_matrix(         scale_m,         sample_step,        size_u/N_u,        size_v/N_v); 
        et_create_translation_matrix(   translation_m,   -0.5 * sample_step * N_samples,    -0.5*size_u,        -0.5*size_v);
    }
    else if (direction==4) {
        dim[1] = N_samples;          dim[2] = N_v;              dim[3] = N_u; 
        et_create_scale_matrix(         scale_m,         sample_step,         size_v/N_v,         size_u/N_u); 
        et_create_translation_matrix(   translation_m,   -0.5 * sample_step * N_samples,      -0.5*size_v,        -0.5*size_u);
    }
    else if (direction==5) {
        dim[1] = N_v;                dim[2] = N_samples;        dim[3] = N_u; 
        et_create_scale_matrix(         scale_m,         size_v/N_v,         sample_step,       size_u/N_u); 
        et_create_translation_matrix(   translation_m,   -0.5*size_v,    -0.5 * sample_step * N_samples,     -0.5*size_u);
    }
    else if (direction==6) {
        dim[1] = N_u;                dim[2] = N_samples;        dim[3] = N_v; 
        et_create_scale_matrix(         scale_m,         size_u/N_u,         sample_step,       size_v/N_v); 
        et_create_translation_matrix(   translation_m,   -0.5*size_u,    -0.5 * sample_step * N_samples,    -0.5*size_v);
    }
    else {
        dim[1] = N_v;                dim[2] = N_u;              dim[3] = N_samples; 
        et_create_scale_matrix(         scale_m,         size_v/N_v,         size_u/N_u,        sample_step); 
        et_create_translation_matrix(   translation_m,   -0.5*size_v,        -0.5*size_u,       -0.5 * sample_step * N_samples);
    }

    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *resampledImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    
    *m = reg_mat44_mul(translation_m, scale_m); 

    //  2) set sform matrix: 
    resampledImage->pixdim[0] = 1;    
    resampledImage->qform_code = 0; 
    resampledImage->sform_code = 1; 
    resampledImage->sto_xyz.m[0][0]=m->m[0][0];     resampledImage->sto_xyz.m[0][1]=m->m[0][1];    resampledImage->sto_xyz.m[0][2]=m->m[0][2];     resampledImage->sto_xyz.m[0][3]=m->m[0][3];
    resampledImage->sto_xyz.m[1][0]=m->m[1][0];     resampledImage->sto_xyz.m[1][1]=m->m[1][1];    resampledImage->sto_xyz.m[1][2]=m->m[1][2];     resampledImage->sto_xyz.m[1][3]=m->m[1][3];
    resampledImage->sto_xyz.m[2][0]=m->m[2][0];     resampledImage->sto_xyz.m[2][1]=m->m[2][1];    resampledImage->sto_xyz.m[2][2]=m->m[2][2];     resampledImage->sto_xyz.m[2][3]=m->m[2][3]; 
    resampledImage->sto_xyz.m[3][0]=m->m[3][0];     resampledImage->sto_xyz.m[3][1]=m->m[3][1];    resampledImage->sto_xyz.m[3][2]=m->m[3][2];     resampledImage->sto_xyz.m[3][3]=m->m[3][3];
    resampledImage->sto_ijk = nifti_mat44_inverse(resampledImage->sto_xyz); 

    //  3) allocate GPU memory
    if(cudaCommon_allocateArrayToDevice<float>(&resampledActivityArray_d, resampledImage->dim)) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough memory on GPU for 'resampledActivityArray_d' - %lu bytes \n",resampledImage->ndim*sizeof(float));
        return niftyrec_error_allocgpu;
    }  
    alloc_record_add(memory_record,(void*)resampledActivityArray_d,ALLOCTYPE_CUDA); 
        
    // -- Attenuation -- 
//    if (attenuationImage != NULL) {
//        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) {
//            alloc_record_destroy(memory_record); 
//            fprintf(stderr,"Not enough memory on GPU for 'attenuationArray_d' - %lu bytes \n",attenuationImage->ndim*sizeof(float));
//            return niftyrec_error_transfergpu;
//        } 
//        alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY); 
//        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) {
//            alloc_record_destroy(memory_record);
//            fprintf(stderr,"Unable to copy %lu bytes of memory to GPU: from 'attenuationImage->data' to 'attenuationArray_d'. \n",attenuationImage->ndim*sizeof(float));
//            return niftyrec_error_transfergpu;
//        }
//    } 

//    // -- Resampled Attenuation -- 
//    if (attenuationImage != NULL) {
//        if(cudaCommon_allocateArrayToDevice<float>(&resampledAttenuationArray_d, resampledImage->dim)) {
//            alloc_record_destroy(memory_record); 
//            fprintf(stderr,"Not enough memory on GPU for 'resampledAttenuationArray_d' - %lu bytes \n",resampledImage->ndim*sizeof(float));
//            return niftyrec_error_allocgpu;
//        }
//        alloc_record_add(memory_record,(void*)resampledAttenuationArray_d,ALLOCTYPE_CUDA);
//    }

    // -- Locations --
    if(cudaMalloc((void **)&locationsArray_d, 3*N_locations*sizeof(unsigned short)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough memory on GPU for 'locationsArray_d' - %lu bytes \n",3*N_locations*sizeof(unsigned short));
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)locationsArray_d,ALLOCTYPE_CUDA);
    
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    gettimeofday(&t1, NULL);
    cudaMemcpy(locationsArray_d, locations, 3*N_locations*sizeof(unsigned short), cudaMemcpyHostToDevice);
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    // -- Projection --  
    gettimeofday(&t1, NULL);
    if(cudaMalloc((void **)&projectionArray_d, N_locations*sizeof(float)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough memory on GPU for 'projectionArray_d' - %lu bytes \n",N_locations*sizeof(float));
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)projectionArray_d,ALLOCTYPE_CUDA);
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 
    gettimeofday(&t1, NULL);
    cudaMemset((void*)projectionArray_d,0.0f,N_locations*sizeof(float));
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    // -- Position field --
    gettimeofday(&t1, NULL);
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, resampledImage->dim)) { 
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough memory on GPU for 'positionFieldImageArray_d' - %lu bytes \n",resampledImage->ndim*sizeof(float));
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA); 
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    // -- Mask -- 
    gettimeofday(&t1, NULL);
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, resampledImage->dim)) {  
         alloc_record_destroy(memory_record); 
         fprintf(stderr,"Not enough memory on GPU for 'mask_d' - %lu bytes \n",resampledImage->ndim*sizeof(int));
         return niftyrec_error_transfergpu;
    } 
    alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);
    int *mask_h=(int *)malloc(resampledImage->nvox*sizeof(int));
    for(int i=0; i<resampledImage->nvox; i++) mask_h[i]=i;
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    

    gettimeofday(&t1, NULL);    
    if (cudaMemcpy(mask_d, mask_h, resampledImage->nvox*sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Unable to transfer %lu bytes to GPU: from 'mask_h' to 'mask_d'. \n",resampledImage->nvox*sizeof(int));
        return niftyrec_error_transfergpu;
    }
    free(mask_h);
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    // -- Transformation matrix -- 
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_HOST);

/*
//DEBUG: 
    dim[0]    = 3; 
    dim[1]    = N_u; 
    dim[2]    = N_v;        
    dim[3]    = N_samples; 
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *debugImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    debugImage->dim[0]=dim[0]; debugImage->dim[1]=dim[1]; debugImage->dim[2]=dim[2]; debugImage->dim[3]=dim[3]; debugImage->dim[4]=dim[4]; debugImage->dim[5]=dim[5]; debugImage->dim[6]=dim[6]; debugImage->dim[7]=dim[7]; 
    debugImage->nx=dim[1]; debugImage->ny=dim[2]; debugImage->nz=dim[3]; debugImage->nt=dim[4];  
    debugImage->nvox = dim[1]*dim[2]*dim[3]*dim[4];
    debugImage->data = (float *)malloc(debugImage->nvox*sizeof(float));  
    debugImage->pixdim[0] = 1;    
    debugImage->qform_code = 0; 
    debugImage->sform_code = 1; 
    debugImage->sto_xyz.m[0][0]=1;     debugImage->sto_xyz.m[0][1]=0;    debugImage->sto_xyz.m[0][2]=0;     debugImage->sto_xyz.m[0][3]=0;
    debugImage->sto_xyz.m[1][0]=0;     debugImage->sto_xyz.m[1][1]=1;    debugImage->sto_xyz.m[1][2]=0;     debugImage->sto_xyz.m[1][3]=0;
    debugImage->sto_xyz.m[2][0]=0;     debugImage->sto_xyz.m[2][1]=0;    debugImage->sto_xyz.m[2][2]=1;     debugImage->sto_xyz.m[2][3]=0; 
    debugImage->sto_xyz.m[3][0]=0;     debugImage->sto_xyz.m[3][1]=0;    debugImage->sto_xyz.m[3][2]=0;     debugImage->sto_xyz.m[3][3]=1;
    debugImage->sto_ijk = nifti_mat44_inverse(debugImage->sto_xyz); 

    dim[0]    = 4; 
    dim[1]    = N_u; 
    dim[2]    = N_v;        
    dim[3]    = N_samples; 
    dim[4]    = 4; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *tfieldImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    tfieldImage->dim[0]=dim[0]; tfieldImage->dim[1]=dim[1]; tfieldImage->dim[2]=dim[2]; tfieldImage->dim[3]=dim[3]; tfieldImage->dim[4]=dim[4]; tfieldImage->dim[5]=dim[5]; tfieldImage->dim[6]=dim[6]; tfieldImage->dim[7]=dim[7]; 
    tfieldImage->nx=dim[1]; tfieldImage->ny=dim[2]; tfieldImage->nz=dim[3]; tfieldImage->nt=dim[4];  
    tfieldImage->nvox = dim[1]*dim[2]*dim[3]*dim[4];
    tfieldImage->data = (float *)malloc(4*tfieldImage->nvox*sizeof(float));  
    tfieldImage->pixdim[0] = 1;    
    tfieldImage->qform_code = 0; 
    tfieldImage->sform_code = 1; 
    tfieldImage->sto_xyz.m[0][0]=1;     tfieldImage->sto_xyz.m[0][1]=0;    tfieldImage->sto_xyz.m[0][2]=0;     tfieldImage->sto_xyz.m[0][3]=0;
    tfieldImage->sto_xyz.m[1][0]=0;     tfieldImage->sto_xyz.m[1][1]=1;    tfieldImage->sto_xyz.m[1][2]=0;     tfieldImage->sto_xyz.m[1][3]=0;
    tfieldImage->sto_xyz.m[2][0]=0;     tfieldImage->sto_xyz.m[2][1]=0;    tfieldImage->sto_xyz.m[2][2]=1;     tfieldImage->sto_xyz.m[2][3]=0; 
    tfieldImage->sto_xyz.m[3][0]=0;     tfieldImage->sto_xyz.m[3][1]=0;    tfieldImage->sto_xyz.m[3][2]=0;     tfieldImage->sto_xyz.m[3][3]=1;
    tfieldImage->sto_ijk = nifti_mat44_inverse(tfieldImage->sto_xyz); 
*/
    
    gettimeofday(&t1, NULL);
    
    /* Compute projections */ 
    for (unsigned int index_axial=0; index_axial<N_axial; index_axial++) { 
       for (unsigned int index_azim=0; index_azim<N_azimuthal; index_azim++) { 

          unsigned int offsets_index = index_azim*N_axial + index_axial;
          unsigned int is_active = active[offsets_index]; 
          unsigned int index_axial_matrix     = index_azim*N_axial + index_axial;
          unsigned int index_azim_matrix = index_azim*N_axial + index_axial + N_axial*N_azimuthal;
          // Determine the offset of locations array for the current detector plane                    
          int offset = offsets[offsets_index]; 
           
//          fprintf(stdout,"Axial: %d  Azimuthal: %d (offset index: %d) ACTIVE: %d \n",angles[index_axial_matrix],angles[index_azim_matrix],offsets_index,is_active);
//fprintf(stderr, "ROTATION:  I_axial: %d   I_azim: %d   axial: %3.2f   azim: (0.5*PI)+%3.2f  active: %d   I_offset: %d   offset: %d \n", index_axial, index_azim, angles[index_axial_matrix], angles[index_azim_matrix], is_active, offsets_index, offset ); 
           
          if (is_active) {
               
            gettimeofday(&t1, NULL); 
   
            // Apply affine //
            if (direction==1)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.5*PI + angles[index_azim_matrix],
                                            0,
                                            0.5*PI + angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            XZY_ROTATION);
            else if (direction==2)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.0,
                                            0.5*PI + angles[index_azim_matrix],
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YZX_ROTATION);
            else if (direction==3)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.0, 
                                            angles[index_azim_matrix], 
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YZX_ROTATION);
            else if (direction==4)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.5*PI, 
                                            angles[index_azim_matrix], 
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            XYZ_ROTATION);
            else if (direction==5)
                et_create_rotation_matrix(  affineTransformation, 
                                            angles[index_azim_matrix],
                                            0.5*PI, 
                                            0.5*PI + angles[index_axial_matrix],   // should it be -0.5*PI+.. 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YXZ_ROTATION);          
            else if (direction==6)
                et_create_rotation_matrix(  affineTransformation, 
                                            angles[index_azim_matrix], 
                                            0.0,
                                            0.5*PI + angles[index_axial_matrix],   // should it be -0.5*PI+.. 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            XZY_ROTATION);
            else                            
                et_create_rotation_matrix(  affineTransformation, 
                                            0.0,
                                            0.5*PI + angles[index_azim_matrix],
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YZX_ROTATION);
    

            reg_affine_positionField_gpu(   affineTransformation, 
                                            resampledImage, 
                                            &positionFieldImageArray_d);
              
            gettimeofday(&t2, NULL); 
            T2_rotate += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 
            
            gettimeofday(&t1, NULL);
            // Resample the activity image // 
            if (activityImage != NULL)
                reg_resampleSourceImage_gpu(resampledImage,
                                            activityImage,
                                            &resampledActivityArray_d,
                                            &activityArray_d,
                                            &positionFieldImageArray_d,
                                            &mask_d,
                                            resampledImage->nvox,
                                            background,
                                            0);

            gettimeofday(&t2, NULL); 
            T3_resample += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 

//            // Resample the attenuation map // 
//            if (attenuationImage != NULL)
//                reg_resampleSourceImage_gpu( 
//                                            resampledImage,
//                                            attenuationImage,
//                                            &resampledAttenuationArray_d,
//                                            &attenuationArray_d,
//                                            &positionFieldImageArray_d,
//                                            &mask_d,
//                                            resampledImage->nvox,
//                                            background_attenuation,
//                                            0);




// DEBUG: 
/*
if (index_axial==0 && index_azim==0) {
	        fprintf(stderr,"=======   0) Activity volume:   ======\n");
	        for (int i=0; i<4; i++)
			    fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",activityImage->sto_xyz.m[i][0],activityImage->sto_xyz.m[i][1],activityImage->sto_xyz.m[i][2],activityImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

	        fprintf(stderr,"=======   0) Resampled volume:  ======\n");
	        for (int i=0; i<4; i++)
			    fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",resampledImage->sto_xyz.m[i][0],resampledImage->sto_xyz.m[i][1],resampledImage->sto_xyz.m[i][2],resampledImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

            fprintf(stderr,"=======   0) Transformation:    ======\n");
            for (int i=0; i<4; i++)
    			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",affineTransformation->m[i][0],affineTransformation->m[i][1],affineTransformation->m[i][2],affineTransformation->m[i][3]);
            fprintf(stderr,"=====================================\n"); 
            
    if (cudaMemcpy(((void *)debugImage->data), resampledActivityArray_d, N_u*N_v*N_samples*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    }

    nifti_set_filenames(debugImage,"/Users/spedemon/Desktop/debug00.nii",0,0); 
    nifti_image_write(debugImage); 

    if (cudaMemcpy(((void *)tfieldImage->data), positionFieldImageArray_d, N_u*N_v*N_samples*sizeof(float)*4, cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    }

    nifti_set_filenames(tfieldImage,"/Users/spedemon/Desktop/tfield00.nii",0,0); 
    nifti_image_write(tfieldImage); 
}
if (index_axial==10 && index_azim==0) {
	        fprintf(stderr,"=======   10) Activity volume:   ======\n");
	        for (int i=0; i<4; i++)
			    fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",activityImage->sto_xyz.m[i][0],activityImage->sto_xyz.m[i][1],activityImage->sto_xyz.m[i][2],activityImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

	        fprintf(stderr,"=======   10) Resampled volume:  ======\n");
	        for (int i=0; i<4; i++)
			    fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",resampledImage->sto_xyz.m[i][0],resampledImage->sto_xyz.m[i][1],resampledImage->sto_xyz.m[i][2],resampledImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

            fprintf(stderr,"=======   10) Transformation:    ======\n");
            for (int i=0; i<4; i++)
    			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",affineTransformation->m[i][0],affineTransformation->m[i][1],affineTransformation->m[i][2],affineTransformation->m[i][3]);
            fprintf(stderr,"=====================================\n"); 
            
    if (cudaMemcpy(((void *)debugImage->data), resampledActivityArray_d, N_u*N_v*N_samples*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    } 
    nifti_set_filenames(debugImage,"/Users/spedemon/Desktop/debug10.nii",0,0); 
    nifti_image_write(debugImage); 
    
    if (cudaMemcpy(((void *)tfieldImage->data), positionFieldImageArray_d, N_u*N_v*N_samples*sizeof(float)*4, cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    }
    nifti_set_filenames(tfieldImage,"/Users/spedemon/Desktop/tfield10.nii",0,0); 
    nifti_image_write(tfieldImage); 
}
*/
            // Apply Depth Dependent Point Spread Function //
            // FIXME: implement

            // Integrate along active lines // 
            
            // 2) determine the number of active locations for the current detector plane 
            
            unsigned int N_locations_plane = 0;
            if ((index_azim==(N_azimuthal-1)) && (index_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                unsigned int offsets_index_next = 0;
                unsigned int index_axial_next   = index_axial; 
                unsigned int index_azim_next    = index_azim; 
                if ( index_azim == (N_azimuthal-1) ) {
                     index_axial_next += 1;
                     index_azim_next   = 0; 
                }
                else
                     index_azim_next   += 1;
                offsets_index_next = index_azim_next*N_axial + index_axial_next;
                N_locations_plane = offsets[offsets_index_next] - offset; 
                //fprintf(stderr," -- X index_axial: %d   index_azimuthal: %d   index: %d    offsets[index]: %d    offset: %d \n", index_axial, index_azim, offsets_index_next, offsets[offsets_index_next], offset);
            } 

            // 3) take the pointer to the array of locations for the current detector plane  
            unsigned short *current_locationsArray_d  = locationsArray_d  + 3*offset;   
            float          *current_projectionArray_d = projectionArray_d + offset;  

            /*
            fprintf(stderr," --5 N_locations_plane:  %d \n",N_locations_plane);
            fprintf(stderr," --5 N_locations:        %d \n",N_locations);
            fprintf(stderr," --5 N_u:                %d \n",N_u);
            fprintf(stderr," --5 N_v:                %d \n",N_v);
            fprintf(stderr," --5 N_samples:          %d \n",N_samples);
            fprintf(stderr," --5 direction:          %d \n",direction);
            fprintf(stderr," --5 block_size:         %d \n",block_size);
            fprintf(stderr," --5 offset:             %d \n",offset);
            fprintf(stderr," --5 locationsArray_d:   %dx%d \n",3*N_locations,sizeof(unsigned short) );
            fprintf(stderr," --5 projectionArray_d:  %dx%d \n",N_locations,sizeof(float) );
            fprintf(stderr," --5 locations at start:     %d %d %d   %d %d %d   %d %d %d   %d %d %d   %d %d %d\n", locations[0],locations[1],locations[2],locations[3],locations[4],locations[5],locations[6],locations[7],locations[8],locations[9],locations[10],locations[11],locations[12],locations[13],locations[14]); 
            fprintf(stderr," --5 locations at offset:    %d %d %d   %d %d %d   %d %d %d   %d %d %d   %d %d %d\n", locations[3*offset],locations[3*offset+1],locations[3*offset+2],locations[3*offset+3],locations[3*offset+4],locations[3*offset+5],locations[3*offset+6],locations[3*offset+7],locations[3*offset+8],locations[3*offset+9],locations[3*offset+10],locations[3*offset+11],locations[3*offset+12],locations[3*offset+13],locations[3*offset+14]); 
            fprintf(stderr," --5 locations at end:       %d %d %d   %d %d %d   %d %d %d\n", locations[3*(N_locations-3)], locations[3*(N_locations-3)+1], locations[3*(N_locations-3)+2], locations[3*(N_locations-3)+3], locations[3*(N_locations-3)+4], locations[3*(N_locations-3)+5],locations[3*(N_locations-3)+6],locations[3*(N_locations-3)+7],locations[3*(N_locations-3)+8]);
            */
            
            gettimeofday(&t1, NULL);
            // 4) launch the kernel that computes the projections 
            pet_line_integral_compressed_gpu(
                                            resampledActivityArray_d, 
                                            resampledAttenuationArray_d, 
                                            current_projectionArray_d, 
                                            current_locationsArray_d, 
                                            N_locations_plane, 
                                            N_u, 
                                            N_v,
                                            N_samples,
                                            direction,
                                            block_size); 
            
            gettimeofday(&t2, NULL); 
            T4_integral += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 

          }
          //else 
          //    fprintf(stderr,"Axial: %d  Azimuthal: %d  NOT ACTIVE \n",index_axial,index_azim);
        }
    }

    gettimeofday(&t1, NULL); 
    /* Transfer result back to host */    
    if (cudaMemcpy(projection, projectionArray_d, N_locations*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;
    }
    
    gettimeofday(&t2, NULL);
    T5_transfer_to_host = 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 

    /* Truncate negative values: small negative values may be found due to FFT and IFFT */
    if (truncate_negative_values) {
        for (int i=0; i<N_locations; i++)
            if (projection[i] < 0)
                projection[i] = 0;
    }

    gettimeofday(&t1, NULL); 
    /*Free*/
    if ( alloc_record_destroy(memory_record) ) 
            return niftyrec_error_alloccpu;   
    gettimeofday(&t2, NULL);
    T6_free = 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    time_profiling[0]=T0_transfer_to_gpu; time_profiling[1]=T1_alloc;    time_profiling[2]=T2_rotate; 
    time_profiling[3]=T3_resample;        time_profiling[4]=T4_integral; time_profiling[5]=T5_transfer_to_host; 
    time_profiling[6]=T6_free; 
    
    return niftyrec_success; 
}
#endif 




#ifdef _USE_CUDA
unsigned int PET_project_compressed_gpu_test(nifti_image *activityImage, nifti_image *attenuationImage, float *projection, 
               int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal)
{
    fprintf(stdout,"Projecting compressed PET (GPU) ...\n"); 

    /* define variables - note that *_d indicates variables that reside on the cuda (d)evice, they are of two types: cudaArray (for texture fetching) or not cudaArray */
    float          *activityArray_d             = NULL;     //stores input activity, makes use of GPU texture fetch unit
    float          *attenuationArray_d          = NULL;     //stores input attenuation coefficients, makes use of GPU texture fetch unit
    unsigned short *locationsArray_d            = NULL;     //stores the locations array (compressed projection data structure)
    float          *projectionArray_d           = NULL;     //stores the projection data (output)
    
    int N_u = activityImage->nx; 
    int N_v = activityImage->ny; 
    int N_samples = activityImage->nz; 
    
    /* Allocate arrays on the device and transfer data to the device */
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);  //memory handling: keep a record of the allocated memory 

    // -- Activity    
    if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    }  
    alloc_record_add(memory_record,(void*)activityArray_d,ALLOCTYPE_CUDA); 
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_transfergpu;
    }
            
    // -- Attenuation -- 
    if (attenuationImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;
        }  
        alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA); 
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_transfergpu;
        }
    } 

    // -- Locations --
    if(cudaMalloc((void **)&locationsArray_d, 3*N_locations*sizeof(unsigned short)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)locationsArray_d,ALLOCTYPE_CUDA);
    cudaMemcpy(locationsArray_d, locations, 3*N_locations*sizeof(unsigned short), cudaMemcpyHostToDevice);

    // -- Projection --    
    if(cudaMalloc((void **)&projectionArray_d, N_locations*sizeof(float)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)projectionArray_d,ALLOCTYPE_CUDA);




    /* Compute projections */ 
    for (unsigned int index_azim=0; index_azim<N_azimuthal; index_azim++) { 
        for (unsigned int index_axial=0; index_axial<N_axial; index_axial++) { 
          if (active[index_azim*N_axial + index_axial]>0) {
          
            // Determine the offset of locations array for the current detector plane 
            int offset = offsets[index_azim*N_axial + index_axial]; 

            // Apply Depth Dependent Point Spread Function //
            // FIXME: implement

            // Integrate along active lines // 
            // 1) determine the number of active locations for the current detector plane 
            unsigned int N_locations_plane = 0;
            if ((index_azim==(N_azimuthal-1)) && (index_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                N_locations_plane = offsets[index_azim*N_axial + index_axial + 1] - offset; 
            } 
            
            // 2) take the pointer to the array of locations for the current detector plane  
            unsigned short *current_locationsArray_d  = locationsArray_d  + 3*offset; 
            float          *current_projectionArray_d = projectionArray_d + offset;   

            // 3) launch the kernel that computes the projections 
            pet_line_integral_compressed_gpu(
                                            activityArray_d, 
                                            attenuationArray_d, 
                                            current_projectionArray_d, 
                                            current_locationsArray_d, 
                                            N_locations_plane, 
                                            N_u, 
                                            N_v,
                                            N_samples,
                                            1,
                                            256 ); 
          }
        }
    }

    /* Transfer result back to host */
    if (cudaMemcpy(projection, projectionArray_d, N_locations*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;
    }

    /*Free*/
    if ( alloc_record_destroy(memory_record) ) 
            return niftyrec_error_alloccpu;     
    return niftyrec_success; 
}
#endif



unsigned int PET_project_compressed_cpu(nifti_image *activityImage, nifti_image *attenuationImage, float *projection, 
               int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
               float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
               unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int truncate_negative_values, 
               unsigned int direction)
{
    fprintf(stdout,"Projecting compressed PET ...\n");
    return 0; 
}



unsigned int next_multiple_of_128(unsigned int N)
{
    return ceil(N/128.0)*128;
}
unsigned int next_multiple_of_64(unsigned int N)
{
    return ceil(N/64.0)*64;
}



#ifdef _USE_CUDA
unsigned int PET_backproject_compressed_gpu(nifti_image *backprojectionImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
                float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
                unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int direction, unsigned int block_size, float *time_profiling)
{
    // performance parameters: default values ('direction' only affects performance, identical results with all 6 values; 'block_size' only affects performance - it's the number of parallel threads per block)
//    if (direction==0)
//        direction  = 2;
    if (block_size==0)
        block_size = 512; 
//fprintf(stderr,"Direction: %d \n",direction);
        
    struct timeval t1, t2; 
    double T0_transfer_to_gpu, T1_alloc, T2_rotate, T3_resample, T4_integral, T5_transfer_to_host, T6_free, T7_accumulate, T8_clear_memory, T9_copy_texture;
    T0_transfer_to_gpu = 0.0; 
    T1_alloc = 0.0; 
    T2_rotate = 0.0; 
    T3_resample = 0.0;
    T4_integral = 0.0; 
    T5_transfer_to_host = 0.0;
    T6_free = 0.0; 
    T7_accumulate = 0.0;
    T8_clear_memory = 0.0;
    T9_copy_texture = 0.0;
    gettimeofday(&t1, NULL);

    
    unsigned int N_u_pitched = next_multiple_of_128(N_u); 
    unsigned int N_v_pitched = next_multiple_of_128(N_v); 

    /* define variables - note that *_d indicates variables that reside on the cuda (d)evice, they are of two types: cudaArray (for texture fetching) or not cudaArray */
    cudaArray      *backprojectionArray_d          = NULL;     //stores the backprojection (output) 
    cudaArray      *attenuationArray_d             = NULL;     //stores input attenuation coefficients, makes use of GPU texture fetch unit
    float          *tempBackprojectionArray_d      = NULL;     //temporarily stores activity aligned with one of the detector plane  
    float          *resampledBackprojectionArray_d = NULL;
    float          *resampledAttenuationArray_d    = NULL;     //temporarily stores attenuation coefficients aligned with one of the detector planes
    float          *attenuationPlaneArray_d        = NULL;     
    float          *integralBackprojectionArray_d  = NULL; 
    float4         *positionFieldImageArray_d      = NULL;     //temporarily stores position field for rotation of activity and attenuation 
    float4         *positionFieldResampledArray_d  = NULL;     //temporarily stores position field for rotation of activity and attenuation 
    int            *maskImage_d                    = NULL;     //binary mask that defines active voxels (typically all active)
    int            *maskResampled_d                = NULL;     //binary mask that defines active voxels (typically all active)
    unsigned short *locationsArray_d               = NULL;     //stores the locations array (compressed projection data structure)
    float          *projectionArray_d              = NULL;     //stores the projection data (input)

    /* Allocate arrays on the GPU device and transfer data to the GPU device */
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);  //memory handling: keep a record of the allocated memory 


    // -- Resampled Image -- 
    // This nifti_image defines the (moving) volume that resamples the (static) imaging volume 

    //  1) define sform matrix: 
    mat44 *scale_m       = (mat44 *)calloc(1,sizeof(mat44));
    mat44 *translation_m = (mat44 *)calloc(1,sizeof(mat44));
    mat44 *m             = (mat44 *)calloc(1,sizeof(mat44));
    alloc_record_add(memory_record,(void*)scale_m,ALLOCTYPE_HOST);
    alloc_record_add(memory_record,(void*)translation_m,ALLOCTYPE_HOST);
    alloc_record_add(memory_record,(void*)m,ALLOCTYPE_HOST);
    
    // 2) define resampled-activity nifti image
    int dim[8]; 
    dim[0]    = 3; 
    if (direction==1) {
        dim[1] = N_u_pitched;                dim[2] = N_v_pitched;              dim[3] = N_samples; 
        et_create_scale_matrix(         scale_m,         size_u/N_u,         size_v/N_v,        sample_step); 
        et_create_translation_matrix(   translation_m,   -0.5*size_u,        -0.5*size_v,       -0.5 * sample_step * N_samples);
    }
    else if (direction==2) {
        dim[1] = N_v;                dim[2] = N_u;              dim[3] = N_samples; 
        et_create_scale_matrix(         scale_m,         size_v/N_v,         size_u/N_u,        sample_step); 
        et_create_translation_matrix(   translation_m,   -0.5*size_v,        -0.5*size_u,       -0.5 * sample_step * N_samples);
    }
    else if (direction==3) {
        dim[1] = N_samples;          dim[2] = N_u;              dim[3] = N_v; 
        et_create_scale_matrix(         scale_m,         sample_step,        size_u/N_u,        size_v/N_v); 
        et_create_translation_matrix(   translation_m,   -0.5 * sample_step * N_samples,    -0.5*size_u,        -0.5*size_v);
    }
    else if (direction==4) {
        dim[1] = N_samples;          dim[2] = N_v;              dim[3] = N_u; 
        et_create_scale_matrix(         scale_m,         sample_step,         size_v/N_v,         size_u/N_u); 
        et_create_translation_matrix(   translation_m,   -0.5 * sample_step * N_samples,      -0.5*size_v,        -0.5*size_u);
    }
    else if (direction==5) {
        dim[1] = N_v;                dim[2] = N_samples;        dim[3] = N_u; 
        et_create_scale_matrix(         scale_m,         size_v/N_v,         sample_step,       size_u/N_u); 
        et_create_translation_matrix(   translation_m,   -0.5*size_v,    -0.5 * sample_step * N_samples,     -0.5*size_u);
    }
    else if (direction==6) {
        dim[1] = N_u;                dim[2] = N_samples;        dim[3] = N_v; 
        et_create_scale_matrix(         scale_m,         size_u/N_u,         sample_step,       size_v/N_v); 
        et_create_translation_matrix(   translation_m,   -0.5*size_u,    -0.5 * sample_step * N_samples,    -0.5*size_v);
    }
    else {
        dim[1] = N_v_pitched;        dim[2] = N_u_pitched;              dim[3] = N_samples; 
        et_create_scale_matrix(         scale_m,         size_v/N_v,         size_u/N_u,        sample_step); 
        et_create_translation_matrix(   translation_m,   -0.5*size_v,        -0.5*size_u,       -0.5 * sample_step * N_samples);
    } 
    
    dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1; 
    nifti_image *resampledImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false); 
    
    *m = reg_mat44_mul(translation_m, scale_m); 

    //  2) set sform matrix: 
    resampledImage->pixdim[0] = 1;    
    resampledImage->qform_code = 0; 
    resampledImage->sform_code = 1; 
    resampledImage->sto_xyz.m[0][0]=m->m[0][0];     resampledImage->sto_xyz.m[0][1]=m->m[0][1];    resampledImage->sto_xyz.m[0][2]=m->m[0][2];     resampledImage->sto_xyz.m[0][3]=m->m[0][3];
    resampledImage->sto_xyz.m[1][0]=m->m[1][0];     resampledImage->sto_xyz.m[1][1]=m->m[1][1];    resampledImage->sto_xyz.m[1][2]=m->m[1][2];     resampledImage->sto_xyz.m[1][3]=m->m[1][3];
    resampledImage->sto_xyz.m[2][0]=m->m[2][0];     resampledImage->sto_xyz.m[2][1]=m->m[2][1];    resampledImage->sto_xyz.m[2][2]=m->m[2][2];     resampledImage->sto_xyz.m[2][3]=m->m[2][3]; 
    resampledImage->sto_xyz.m[3][0]=m->m[3][0];     resampledImage->sto_xyz.m[3][1]=m->m[3][1];    resampledImage->sto_xyz.m[3][2]=m->m[3][2];     resampledImage->sto_xyz.m[3][3]=m->m[3][3];
    resampledImage->sto_ijk = nifti_mat44_inverse(resampledImage->sto_xyz); 

    
    // -- Temporary Backprojection --
    // This volume contains the backprojection from one detector plane, in the space of the detector plane; 
    // The content of this volume is copied into backprojectionArray_d, which is bound to texture. Though this requires more memory, the 
    // volume in texture-bound memory can be resampled tri-linearly with hardware interpolation. 
    cudaPitchedPtr temp_backprojection_pitched; 
    cudaExtent     temp_backprojection_extent = make_cudaExtent(sizeof(float)*resampledImage->nx,resampledImage->ny,resampledImage->nz);     
    cudaError_t    cuda_error = cudaMalloc3D(&temp_backprojection_pitched, temp_backprojection_extent);     
    if(cuda_error != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough GPU memory to allocate 'temp_backprojection_pitched'. \n");
        return niftyrec_error_allocgpu;
    }
    tempBackprojectionArray_d = (float*) temp_backprojection_pitched.ptr;
    alloc_record_add(memory_record,(void*)tempBackprojectionArray_d,ALLOCTYPE_CUDA);


    // -- Backprojection -- 
    // This volume contains the backprojection from one detector plane, in the space of the detector plane
    cudaChannelFormatDesc backprojectionArray_d_chdesc = cudaCreateChannelDesc<float>();
    cudaExtent backprojectionArray_d_extent;
    backprojectionArray_d_extent.width  = resampledImage->nx;
    backprojectionArray_d_extent.height = resampledImage->ny;
    backprojectionArray_d_extent.depth  = resampledImage->nz;
    if (cudaMalloc3DArray(&backprojectionArray_d, &backprojectionArray_d_chdesc, backprojectionArray_d_extent) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough GPU memory to allocate 'backprojectionArray_d'. \n");
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)backprojectionArray_d,ALLOCTYPE_CUDA_ARRAY); 


    // -- Resampled Backprojection -- 
    // This volume contains the backprojection for one detector plane, resampled in the space of the imaging volume. 
    if(cudaCommon_allocateArrayToDevice<float>(&resampledBackprojectionArray_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough GPU memory to allocate 'resampledBackprojectionArray_d'. \n");
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)resampledBackprojectionArray_d,ALLOCTYPE_CUDA); 


    // -- Integral Backprojection -- 
    // This volume contains the result: the sum of resampledBackprojectionArray_d for all detector locations. 
    if(cudaCommon_allocateArrayToDevice<float>(&integralBackprojectionArray_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough GPU memory to allocate 'integralBackprojectionArray_d'. \n");
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)integralBackprojectionArray_d,ALLOCTYPE_CUDA); 
    cudaMemset((void*)integralBackprojectionArray_d,0.0f,backprojectionImage->nvox*sizeof(float));

    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);

    // -- Attenuation -- 
    // This volume contains the attenuation map in the space of the imaging volume 
    if (attenuationImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            fprintf(stderr,"Not enough GPU memory to allocate 'attenuationArray_d'. \n");
            return niftyrec_error_allocgpu;
        }
        gettimeofday(&t2, NULL);
        T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
        gettimeofday(&t1, NULL);
        
        alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            fprintf(stderr,"Unable to transfer 'attenuationArray_d' to GPU. \n");
            return niftyrec_error_transfergpu;
        }
        gettimeofday(&t2, NULL);
        T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
        gettimeofday(&t1, NULL);
    }

    // -- Resampled Attenuation -- 
    // Attenuation resampled in the space of a detector plane
    if (attenuationImage != NULL) {
        if(cudaCommon_allocateArrayToDevice<float>(&resampledAttenuationArray_d, resampledImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            fprintf(stderr,"Not enough GPU memory to allocate 'resampledAttenuationArray_d'. \n");
            return niftyrec_error_allocgpu;
        }
        alloc_record_add(memory_record,(void*)resampledAttenuationArray_d,ALLOCTYPE_CUDA);

        gettimeofday(&t2, NULL);
        T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
        gettimeofday(&t1, NULL);
    }
    

    // -- Locations --
    if(cudaMalloc((void **)&locationsArray_d, 3*N_locations*sizeof(unsigned short)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough GPU memory to allocate 'locationsArray_d'. \n");
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)locationsArray_d,ALLOCTYPE_CUDA);
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);
    cudaMemcpy(locationsArray_d, locations, 3*N_locations*sizeof(unsigned short), cudaMemcpyHostToDevice);
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);

    // -- Projection data --    
    if(cudaMalloc((void **)&projectionArray_d, N_locations*sizeof(float)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        fprintf(stderr,"Not enough GPU memory to allocate 'projectionArray_d'. \n");
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)projectionArray_d,ALLOCTYPE_CUDA);
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);
    cudaMemcpy(projectionArray_d, projection, N_locations*sizeof(float), cudaMemcpyHostToDevice);
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);

    
    // -- Position Field Image -- 
    // Position field for resampling from the space of resampledImage (source) to the imaging volume space (target) 
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough GPU memory to allocate 'positionFieldImageArray_d'. \n");
        return niftyrec_error_allocgpu; 
    }
    alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);
    
/*
    // -- Position Field Resampled -- 
    // Position field for resampling from the imaging volume space (source) to the space of resampledImage (target) 
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldResampledArray_d, resampledImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu; 
    }
    alloc_record_add(memory_record,(void*)positionFieldResampledArray_d,ALLOCTYPE_CUDA);
*/ 
   
    // -- Mask Image -- 
    // Mask for resampling: when resampling into imaging volume space, only locations of the target image corresponding to non-zero maskImage_d are sampled 
    if(cudaCommon_allocateArrayToDevice<int>(&maskImage_d, backprojectionImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough GPU memory to allocate 'maskImage_d'. \n");
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)maskImage_d,ALLOCTYPE_CUDA);
    int *maskImage_h=(int *)malloc(backprojectionImage->nvox*sizeof(int)); 
    if (maskImage_h==NULL) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Not enough CPU memory to allocate 'maskImage_h'. \n");
        return niftyrec_error_alloccpu;
    }
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);
    
    for(int i=0; i<backprojectionImage->nvox; i++) maskImage_h[i]=i;
    if (cudaMemcpy(maskImage_d, maskImage_h, backprojectionImage->nvox*sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        fprintf(stderr,"Unable to copy memory to GPU: from '' to 'maskImage_d'. \n");
        return niftyrec_error_transfergpu;
    }
    free(maskImage_h);
    gettimeofday(&t2, NULL);
    T0_transfer_to_gpu += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);

/*
    // -- Mask Resampled -- 
    // Mask for resampling: when resampling into the space of resampledImage, only locations of the target image corresponding to non-zero maskResampled_d are sampled 
    if(cudaCommon_allocateArrayToDevice<int>(&maskResampled_d, resampledImage->dim) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    alloc_record_add(memory_record,(void*)maskResampled_d,ALLOCTYPE_CUDA);
    int *maskResampled_h=(int *)malloc(resampledImage->nvox*sizeof(int)); 
    if (maskResampled_h==NULL) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_alloccpu;
    }
    for(int i=0; i<resampledImage->nvox; i++) maskResampled_h[i]=i;
    if (cudaMemcpy(maskResampled_d, maskResampled_h, resampledImage->nvox*sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_transfergpu;
    }
    free(maskResampled_h);
*/

    // -- Transformation matrix -- 
    // 4X4 transformation matrix: defines the location of a detector plane 
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    if (affineTransformation==NULL) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_alloccpu;
    }
    alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_HOST);
    gettimeofday(&t2, NULL);
    T1_alloc += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    gettimeofday(&t1, NULL);

    // -- Cuda MemCpy parameters -- 
            cudaError_t cuda_status; 
            cudaExtent volumeSize = make_cudaExtent(resampledImage->nx, resampledImage->ny, resampledImage->nz);
            cudaMemcpy3DParms copyparms={0};
            copyparms.extent = volumeSize;
            copyparms.dstArray = backprojectionArray_d;
            copyparms.kind = cudaMemcpyDeviceToDevice; 
            copyparms.srcPtr = temp_backprojection_pitched;




/*
//DEBUG: 
    nifti_image *debugImage = nifti_copy_nim_info(resampledImage); 
    debugImage->data = (float *)malloc(debugImage->nvox*sizeof(float));  
    debugImage->pixdim[0] = 1;    
    debugImage->qform_code = 0; 
    debugImage->sform_code = 1; 
    debugImage->sto_xyz.m[0][0]=1;     debugImage->sto_xyz.m[0][1]=0;    debugImage->sto_xyz.m[0][2]=0;     debugImage->sto_xyz.m[0][3]=0;
    debugImage->sto_xyz.m[1][0]=0;     debugImage->sto_xyz.m[1][1]=1;    debugImage->sto_xyz.m[1][2]=0;     debugImage->sto_xyz.m[1][3]=0;
    debugImage->sto_xyz.m[2][0]=0;     debugImage->sto_xyz.m[2][1]=0;    debugImage->sto_xyz.m[2][2]=1;     debugImage->sto_xyz.m[2][3]=0; 
    debugImage->sto_xyz.m[3][0]=0;     debugImage->sto_xyz.m[3][1]=0;    debugImage->sto_xyz.m[3][2]=0;     debugImage->sto_xyz.m[3][3]=1;
    debugImage->sto_ijk = nifti_mat44_inverse(debugImage->sto_xyz); 

    mat44 *affineTransformation_debug = (mat44 *)calloc(1,sizeof(mat44));
    if (affineTransformation==NULL) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_alloccpu;
    }
    alloc_record_add(memory_record,(void*)affineTransformation_debug,ALLOCTYPE_HOST);
    
    cudaPitchedPtr debug_pitched; 
    cudaExtent     debug_extent = make_cudaExtent(sizeof(float)*resampledImage->nx,resampledImage->ny,resampledImage->nz);     
    cuda_error = cudaMalloc3D(&debug_pitched, debug_extent);     
    if(cuda_error != cudaSuccess) {
        alloc_record_destroy(memory_record);
        return niftyrec_error_allocgpu;
    }
    float *debug_d = (float*) debug_pitched.ptr;
    alloc_record_add(memory_record,(void*)debug_d,ALLOCTYPE_CUDA);

    cudaMemcpy3DParms copyparms1={0};
    copyparms1.extent   = make_cudaExtent(resampledImage->nx, resampledImage->ny, resampledImage->nz);
    copyparms1.srcArray = backprojectionArray_d;
    //copyparms1.srcPtr = temp_backprojection_pitched;
    copyparms1.kind     = cudaMemcpyDeviceToDevice; 
    copyparms1.dstPtr   = debug_pitched;
*/

    gettimeofday(&t1, NULL);
    
    /* Compute backprojection */
    for (unsigned int index_axial=0; index_axial<N_axial; index_axial++) { 
       for (unsigned int index_azim=0; index_azim<N_azimuthal; index_azim++) { 
          unsigned int offsets_index = index_azim*N_axial + index_axial;
          if (active[index_azim*N_axial + index_axial]>0) {
          
            // Determine the offset of locations array for the current detector plane 
            int offset = offsets[index_azim*N_axial + index_axial]; 
            
            unsigned int index_axial_matrix     = index_azim * N_axial + index_axial;
            unsigned int index_azim_matrix = index_azim * N_axial + index_axial + N_axial*N_azimuthal;
              
            // Construct transformation matrix for the current detector plane 
            if (direction==1)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.5*PI + angles[index_azim_matrix],
                                            0,
                                            0.5*PI + angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            XZY_ROTATION);
            else if (direction==2)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.0,
                                            0.5*PI + angles[index_azim_matrix],
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YZX_ROTATION);
            else if (direction==3)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.0, 
                                            angles[index_azim_matrix], 
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YZX_ROTATION);
            else if (direction==4)
                et_create_rotation_matrix(  affineTransformation, 
                                            0.5*PI, 
                                            angles[index_azim_matrix], 
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            XYZ_ROTATION);
            else if (direction==5)
                et_create_rotation_matrix(  affineTransformation, 
                                            angles[index_azim_matrix],
                                            0.5*PI, 
                                            0.5*PI + angles[index_axial_matrix],   // should it be -0.5*PI+.. 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YXZ_ROTATION);          
            else if (direction==6)
                et_create_rotation_matrix(  affineTransformation, 
                                            angles[index_azim_matrix], 
                                            0.0,
                                            0.5*PI + angles[index_axial_matrix],   // should it be -0.5*PI+.. 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            XZY_ROTATION);
            else 
                et_create_rotation_matrix(  affineTransformation, 
                                            0.0,
                                            0.5*PI + angles[index_azim_matrix],
                                            angles[index_axial_matrix], 
                                            0.0, 
                                            0.0, 
                                            0.0, 
                                            YZX_ROTATION);          

            mat44 affineTransformation_inv = nifti_mat44_inverse(*affineTransformation); 
            
            // Transform attenuation map //
            if (attenuationImage != NULL) {
                //no need; could be useful 
            }

            // Ray casting //
            // 1) clear volume
            gettimeofday(&t1, NULL);
            cudaMemset((void*)tempBackprojectionArray_d,0,resampledImage->nvox*sizeof(float));
            gettimeofday(&t2, NULL);
            T8_clear_memory += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
              
            // 2) determine the number of active locations for the current detector plane 
            unsigned int N_locations_plane = 0;
            if ((index_azim==(N_azimuthal-1)) && (index_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                unsigned int offsets_index_next = 0;
                unsigned int index_axial_next   = index_axial; 
                unsigned int index_azim_next    = index_azim; 
                if ( index_azim == (N_azimuthal-1) ) {
                     index_axial_next += 1;
                     index_azim_next   = 0; 
                }
                else
                     index_azim_next   += 1;
                offsets_index_next = index_azim_next*N_axial + index_axial_next;
                N_locations_plane = offsets[offsets_index_next] - offset; 
                //fprintf(stderr," -- X index_axial: %d   index_azimuthal: %d   index: %d    offsets[index]: %d    offset: %d \n", index_axial, index_azim, offsets_index_next, offsets[offsets_index_next], offset);
            } 
            
            // 3) take the pointer to the array of locations for the current detector plane  
            unsigned short *current_locationsArray_d  = locationsArray_d  + 3*offset; 
            float          *current_projectionArray_d = projectionArray_d + offset;   

            gettimeofday(&t1, NULL);
              
            // 4) launch the ray-casting kernel  
            pet_line_backproject_compressed_gpu(
                                            tempBackprojectionArray_d, 
                                            NULL,  // attenuation (pre-integrate instead? ) 
                                            current_projectionArray_d, 
                                            current_locationsArray_d, 
                                            N_locations_plane, 
                                            N_u_pitched, 
                                            N_v_pitched,
                                            N_samples,
                                            direction,
                                            block_size); 

            gettimeofday(&t2, NULL);
            T4_integral += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
              
            // Apply Depth Dependent Point Spread Function //
            // ..

            // Copy the back-projection to texture bound memory (for efficient resampling) //
            gettimeofday(&t1, NULL);
            cuda_status = cudaMemcpy3D(&copyparms);
            if (cuda_status != cudaSuccess) {
                fprintf(stderr, "Error copying to texture-bound memory: %s\n",cudaGetErrorString(cuda_status));
                return 1;
            }
            gettimeofday(&t2, NULL);
            T9_copy_texture += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
              
            // Rotate backprojection //
//cuda_status      = cudaMemcpy3D(&copyparms1);
//if (cuda_status != cudaSuccess) {
//    fprintf(stderr, "##### Error copying from texture-bound memory: %s ##### \n",cudaGetErrorString(cuda_status));
//    return 1;
//}
//et_create_rotation_matrix(  affineTransformation_debug, 
//                                            0.0,
//                                            0.0,
//                                            0.0,
//                                            0.0, 
//                                            0.0, 
//                                            0.0, 
//                                            XZY_ROTATION);
//
            
            gettimeofday(&t1, NULL);
              
            reg_affine_positionField_gpu(     &affineTransformation_inv,
                                              backprojectionImage,
                                              &positionFieldImageArray_d);
             
            gettimeofday(&t2, NULL);
            T2_rotate += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
            gettimeofday(&t1, NULL);
              
            reg_resampleSourceImage_gpu(      backprojectionImage,
                                              resampledImage,
                                              &resampledBackprojectionArray_d,
                                              &backprojectionArray_d,
                                              &positionFieldImageArray_d,
                                              &maskImage_d,
                                              backprojectionImage->nvox,
                                              background,
                                              0);

            gettimeofday(&t2, NULL);
            T3_resample += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
            gettimeofday(&t1, NULL);
              
            // Compute integral backprojection //
            et_accumulate_gpu(                &resampledBackprojectionArray_d,
                                              &integralBackprojectionArray_d,
                                              backprojectionImage );

            gettimeofday(&t2, NULL);
            T7_accumulate += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
            gettimeofday(&t1, NULL);
              
/*
// DEBUG: 
if (index_axial==0 && index_azim==0) {
	        fprintf(stderr,"====  0) Backprojection volume:  =====\n");
	        for (int i=0; i<4; i++)
			    fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",backprojectionImage->sto_xyz.m[i][0],backprojectionImage->sto_xyz.m[i][1],backprojectionImage->sto_xyz.m[i][2],backprojectionImage->sto_xyz.m[i][3]);
            fprintf(stderr,"======================================\n"); 

	        fprintf(stderr,"=======   0) Resampled volume:  ======\n");
	        for (int i=0; i<4; i++)
			    fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",resampledImage->sto_xyz.m[i][0],resampledImage->sto_xyz.m[i][1],resampledImage->sto_xyz.m[i][2],resampledImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

            fprintf(stderr,"=======   0) Transformation:    =====\n");
            for (int i=0; i<4; i++)
    			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",affineTransformation->m[i][0],affineTransformation->m[i][1],affineTransformation->m[i][2],affineTransformation->m[i][3]);
            fprintf(stderr,"=====================================\n"); 

            fprintf(stderr,"=======   0) Transformation_inv:  ===\n");
            for (int i=0; i<4; i++)
    			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",(&affineTransformation_inv)->m[i][0],(&affineTransformation_inv)->m[i][1],(&affineTransformation_inv)->m[i][2],(&affineTransformation_inv)->m[i][3]);
            fprintf(stderr,"=====================================\n"); 
            
    if (cudaMemcpy(((void *)debugImage->data), tempBackprojectionArray_d, N_u*N_v*N_samples*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    }

    nifti_set_filenames(debugImage,"/Users/spedemon/Desktop/debug00.nii",0,0); 
    nifti_image_write(debugImage); 
}
if (index_axial==10 && index_azim==0) {
	        fprintf(stderr,"====  10) Backprojection volume:  ===\n");
	        for (int i=0; i<4; i++)
                fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",backprojectionImage->sto_xyz.m[i][0],backprojectionImage->sto_xyz.m[i][1],backprojectionImage->sto_xyz.m[i][2],backprojectionImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

	        fprintf(stderr,"=======   10) Resampled volume:  ====\n");
	        for (int i=0; i<4; i++)
                fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",resampledImage->sto_xyz.m[i][0],resampledImage->sto_xyz.m[i][1],resampledImage->sto_xyz.m[i][2],resampledImage->sto_xyz.m[i][3]);
            fprintf(stderr,"=====================================\n"); 

            fprintf(stderr,"=======   10) Transformation:    ====\n");
            for (int i=0; i<4; i++)
                fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",affineTransformation->m[i][0],affineTransformation->m[i][1],affineTransformation->m[i][2],affineTransformation->m[i][3]);
            fprintf(stderr,"=====================================\n"); 
          
            fprintf(stderr,"=======   0) Transformation_inv:  ===\n");
            for (int i=0; i<4; i++)
    			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",(&affineTransformation_inv)->m[i][0],(&affineTransformation_inv)->m[i][1],(&affineTransformation_inv)->m[i][2],(&affineTransformation_inv)->m[i][3]);
            fprintf(stderr,"=====================================\n"); 
  
    if (cudaMemcpy(((void *)debugImage->data), tempBackprojectionArray_d, N_u*N_v*N_samples*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    } 
    nifti_set_filenames(debugImage,"/Users/spedemon/Desktop/debug10.nii",0,0); 
    nifti_image_write(debugImage); 
}
*/

          }
        }
    }
    
    gettimeofday(&t1, NULL);
    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(backprojectionImage, &integralBackprojectionArray_d)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu; 
    }
    gettimeofday(&t2, NULL);
    T5_transfer_to_host += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec);
    
    /* Truncate negative values: small negative values may be found due to FFT and IFFT */
    // ..

    // FIXME: this is a quick fix, the back-projection results mirrored, it should not 
//    float *tmp_data = (float *)calloc(backprojectionImage->nvox,sizeof(float)); 
//    alloc_record_add(memory_record,(void*)tmp_data,ALLOCTYPE_HOST); 
//    if (cudaMemcpy(((void *)tmp_data), integralBackprojectionArray_d, backprojectionImage->nvox*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
//        alloc_record_destroy(memory_record); 
//        return niftyrec_error_transfergpu;
//    }
//    float *bdata = (float*) backprojectionImage->data; 
//    for (int z=0; z<backprojectionImage->nz; z++) {
//        for (int x=0; x<backprojectionImage->nx; x++) {
 //           for (int y=0; y<backprojectionImage->ny; y++) {
//                bdata[x + y*backprojectionImage->nx + z*backprojectionImage->nx*backprojectionImage->ny] = tmp_data[x + y*backprojectionImage->nx + (backprojectionImage->nz-z-1)*backprojectionImage->nx*backprojectionImage->ny]; 
//            }
//        }
//    }
        
    

    /*Free*/
    gettimeofday(&t1, NULL);
    if ( alloc_record_destroy(memory_record) ) 
            return niftyrec_error_alloccpu;  
    gettimeofday(&t2, NULL);
    T6_free += 1000000 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec); 
            
    time_profiling[0]=T0_transfer_to_gpu; time_profiling[1]=T1_alloc;    time_profiling[2]=T2_rotate; 
    time_profiling[3]=T3_resample;        time_profiling[4]=T4_integral; time_profiling[5]=T5_transfer_to_host; 
    time_profiling[6]=T6_free;            time_profiling[7]=T7_accumulate; time_profiling[8]=T8_clear_memory;
    time_profiling[9]=T9_copy_texture;           
    
    return niftyrec_success; 
}
#endif 













unsigned int PET_backproject_compressed_cpu(nifti_image *backprojectionImage, nifti_image *attenuationImage, float *projection, 
                int *offsets, unsigned short *locations, unsigned int *active, unsigned int N_locations, unsigned int N_axial, unsigned int N_azimuthal, 
                float *angles, unsigned int N_u, unsigned int N_v, float size_u, float size_v, 
                unsigned int N_samples, float sample_step, float background, float background_attenuation, unsigned int direction) 
{
    return 0; 
}







unsigned int et_spherical_phantom(nifti_image *phantomImage, float centerx, float centery, float centerz, float radius, float inner_value, float outer_value)
{
    float *image = (float*) phantomImage->data; 
    float dx = phantomImage->pixdim[1]; 
    float dy = phantomImage->pixdim[2]; 
    float dz = phantomImage->pixdim[3]; 
    float distance_sq;
    float radius_sq = radius*radius; 

    for (unsigned int x=0; x<phantomImage->nx; x++) {
        for (unsigned int y=0; y<phantomImage->ny; y++) {
            for (unsigned int z=0; z<phantomImage->nz; z++) { 
                distance_sq = (x*dx - centerx)*(x*dx - centerx) + (y*dy - centery)*(y*dy - centery) + (z*dz - centerz)*(z*dz - centerz); 
                if (distance_sq<=radius_sq) {
                    image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = inner_value; 
                }
                else {
                    image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = outer_value; 
                }
                
            }
        }
    }
    return 0; 
}


unsigned int et_cylindrical_phantom(nifti_image *phantomImage, float centerx, float centery, float centerz, float radius, float length, unsigned int axis, float inner_value, float outer_value)
{
    float *image = (float*) phantomImage->data; 
    float dx = phantomImage->pixdim[1]; 
    float dy = phantomImage->pixdim[2]; 
    float dz = phantomImage->pixdim[3]; 
    float distance_sq;
    float radius_sq = radius*radius; 

    // if around X axis: 
    if (axis==0) { 
        for (unsigned int x=0; x<phantomImage->nx; x++) {
            for (unsigned int y=0; y<phantomImage->ny; y++) {
                for (unsigned int z=0; z<phantomImage->nz; z++) { 
                    if ((x*dx - centerx)*(x*dx - centerx) <= (0.5*length)*(0.5*length)) {
                        distance_sq = (y*dy - centery)*(y*dy - centery) + (z*dz - centerz)*(z*dz - centerz); 
                        if (distance_sq<=radius_sq) {
                            image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = inner_value; 
                        }
                        else {
                            image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = outer_value; 
                        }
                    }
                }
            }
        }
    }
    // Y 
    else if (axis==1) { 
        for (unsigned int x=0; x<phantomImage->nx; x++) {
            for (unsigned int y=0; y<phantomImage->ny; y++) {
                for (unsigned int z=0; z<phantomImage->nz; z++) { 
                    if ((y*dy - centery)*(y*dy - centery) <= (0.5*length)*(0.5*length)) {
                        distance_sq = (x*dx - centerx)*(x*dx - centerx) + (z*dz - centerz)*(z*dz - centerz); 
                        if (distance_sq<=radius_sq) {
                            image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = inner_value; 
                        }
                        else {
                            image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = outer_value; 
                        }
                    }
                }
            }
        }
    }
    // Z
    else if (axis==2) { 
        for (unsigned int x=0; x<phantomImage->nx; x++) {
            for (unsigned int y=0; y<phantomImage->ny; y++) {
                for (unsigned int z=0; z<phantomImage->nz; z++) { 
                    if ((z*dz - centerz)*(z*dz - centerz) <= (0.5*length)*(0.5*length)) {
                        distance_sq = (x*dx - centerx)*(x*dx - centerx) + (y*dy - centery)*(y*dy - centery); 
                        if (distance_sq<=radius_sq) {
                            image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = inner_value; 
                        }
                        else {
                            image[x + y*phantomImage->nx + z*phantomImage->nx*phantomImage->ny] = outer_value; 
                        }
                    }
                }
            }
        }
    }
    return 0; 
}


unsigned int et_spheres_ring_phantom(nifti_image *phantomImage, float centerx, float centery, float centerz, float ring_radius, float min_sphere_radius, float max_sphere_radius, unsigned int N_spheres, float inner_value, float outer_value, float taper, unsigned int ring_axis)
{    
    if (ring_axis > 2) {
        return 1; 
    }
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);
    
    // 1) allocate temporary image of the same size as the output image
    nifti_image *tmpImage = nifti_copy_nim_info(phantomImage);
    tmpImage->data = (void*) calloc(tmpImage->nvox, tmpImage->nbyper); 
    alloc_record_add(memory_record,tmpImage,ALLOCTYPE_NIFTI); 

    if (inner_value == 0.0) 
        inner_value = inner_value + 1e-12; 
    float *phantom_data = (float*) phantomImage->data;
    float *tmp_data = (float*) tmpImage->data;
    float c[4]; float c_rotated[4];
    mat44 *rotation_m = (mat44 *)calloc(1,sizeof(mat44)); 
    alloc_record_add(memory_record,(void*)rotation_m,ALLOCTYPE_HOST); 
    float angle; 
    float current_taper; 
     
    // 2) make one sphere at the time and sum into the output image 
    for (int i=0; i<N_spheres; i++) { 
        // Compute the center of the sphere:
        angle         = i*2*PI/ (float) N_spheres; 
        current_taper = -0.5*taper + i*taper/(float)(N_spheres-1);

        if (ring_axis == 0) {
            // Ring around the X axis: 
            c[0] = centerx + current_taper; 
            c[1] = centery + ring_radius; 
            c[2] = centerz; 
            c[3] = 1;
            et_create_rotation_matrix(rotation_m,   angle, 0, 0,   centerx, centery, centerz, XYZ_ROTATION); 
            reg_mat44_mul(rotation_m, c, c_rotated); 
        }
        else if (ring_axis == 1) {
            // Ring around the Y axis: 
            c[0] = centerx; 
            c[1] = centery + current_taper; 
            c[2] = centerz + ring_radius; 
            c[3] = 1;
            et_create_rotation_matrix(rotation_m,   0, angle, 0,   centerx, centery, centerz, XYZ_ROTATION); 
            reg_mat44_mul(rotation_m, c, c_rotated); 
        }
        else if (ring_axis == 2) { 
            // Ring around the Z axis: 
            c[0] = centerx + ring_radius; 
            c[1] = centery; 
            c[2] = centerz + current_taper; 
            c[3] = 1;
            et_create_rotation_matrix(rotation_m,   0, 0, angle,   centerx, centery, centerz, XYZ_ROTATION); 
            reg_mat44_mul(rotation_m, c, c_rotated); 
        }
        // The radius of the sphere increases linearly from min_sphere_radius to max_sphere_radius:
        float r  = min_sphere_radius + i*(max_sphere_radius-min_sphere_radius)/(N_spheres-1); 
        
        // Make sphere: 
        et_spherical_phantom(tmpImage, c_rotated[0], c_rotated[1], c_rotated[2], r, inner_value, 0.0); 
        
        // Set value inside of the spheres (inner_value): 
        for (int j=0; j<phantomImage->nvox; j++) {
            // This sets the value to 'inner_value' also in the regions where two or more spheres eventually intersect: 
            phantom_data[j] += (1-bool(phantom_data[j])) * tmp_data[j];  
        }
        memset((void*) tmp_data, 0, sizeof(float)*tmpImage->nvox); 
    }

    // 3) Set outer value in all voxels that are set to zero (note that if inner_value is set to 0, EPS is added)
    for (int i=0; i<phantomImage->nvox; i++) {
        if (phantom_data[i]==0.0) 
            phantom_data[i] = outer_value;
    }

    /*Free*/
    if ( alloc_record_destroy(memory_record) ) 
            return niftyrec_error_alloccpu;     
    return niftyrec_success; 
}




unsigned int ET_box_to_grid(nifti_image *image, mat44 *affine)
{
    unsigned int status = niftyrec_success; 
    reg_affine_positionField(       affine,
                                    image,
                                    image );
    return status; 
}



unsigned int tr_resample_grid_cpu(nifti_image *resampled, nifti_image *image, nifti_image *grid, float background, unsigned int interpolation_mode)
{
    unsigned int status = niftyrec_success; 
    return status; 
}

#ifdef _USE_CUDA
unsigned int tr_resample_grid_gpu(nifti_image *resampled, nifti_image *image, nifti_image *grid, float background, unsigned int interpolation_mode)
{
    /* initialise the cuda arrays */
    cudaArray *imageArray_d;
    float     *resampledArray_d;
    float4    *gridArray_d;
    int       *mask_d; 

    /* Allocate arrays on the device */
    if(cudaCommon_allocateArrayToDevice<float>(&imageArray_d, image->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&resampledArray_d, resampled->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float4>(&gridArray_d, resampled->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<int>(&mask_d, resampled->dim)) return 1;

    /* Transfer data from the host to the device */
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&imageArray_d,image)) return 1;
    if(cudaCommon_transferNiftiToArrayOnDevice<float4>(&gridArray_d,grid)) return 1;
    int *mask_h=(int *)malloc(resampled->nvox*sizeof(int));
    for(int i=0; i<resampled->nvox; i++) mask_h[i]=i;
    cudaMemcpy(mask_d, mask_h, resampled->nvox*sizeof(int), cudaMemcpyHostToDevice);
    free(mask_h);

   /* Resample the source image */
    reg_resampleSourceImage_gpu(    
                    resampled,
                    image,
                    &resampledArray_d,
                    &imageArray_d,
                    &gridArray_d,
                    &mask_d,
                    resampled->nvox,
                    background, 
                    interpolation_mode);
    /* Transfer result back to host */
    if(cudaCommon_transferFromDeviceToNifti(resampled, &resampledArray_d)) return 1;
    /*Free*/
    cudaCommon_free(&imageArray_d);
    cudaCommon_free((void **)&resampledArray_d);
    cudaCommon_free((void **)&mask_d);
    cudaCommon_free((void **)&gridArray_d);

    return niftyrec_success;
}
#endif




#ifdef _USE_CUDA
unsigned int tr_transform_grid_gpu(nifti_image *transformed_grid, nifti_image *grid, mat44 *affine)
{

    return niftyrec_success;
}
#endif


unsigned int tr_transform_grid_cpu(nifti_image *transformed_grid, nifti_image *grid, mat44 *affine)
{
    return niftyrec_success;
    
}



/*
  Compress and uncompress projection data 
*/

//Uncompress a static projection
unsigned int pet_uncompress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection) 
{
    int status=STATUS_SUCCESS; 

    unsigned int N_locations = *n_locations; 
    unsigned int N_axial     = *n_axial; 
    unsigned int N_azimuthal = * n_azimuthal; 
    unsigned int N_u = *n_u; 
    unsigned int N_v = *n_v; 
    
    unsigned int bin_axial; 
    unsigned int bin_azimuthal; 
    unsigned int bin_u, bin_v, bin_w, i; 
    unsigned int offset; 
    unsigned int N_locations_plane; 
    
    unsigned int offsets_index; 
    unsigned int offsets_index_next; 
    unsigned int bin_azimuthal_next; 
    unsigned int bin_axial_next; 
    
    // for each detector plane 
    for (bin_axial=0; bin_axial < N_axial; bin_axial++) { 
        for (bin_azimuthal=0; bin_azimuthal < N_azimuthal; bin_azimuthal++) { 
            offsets_index = bin_azimuthal*N_axial + bin_axial;
            offset = offsets_matrix[offsets_index]; 
            // the following decompresses the data for a single detector plane: 
            // 1) determine the offset of the data and the number of active locations 
            //fprintf(stdout,"Uncompressing:  axial: %d   azimuthal: %d    offset: %d \n",bin_axial, bin_azimuthal, offset); 
            
            if ((bin_azimuthal==(N_azimuthal-1)) && (bin_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                offsets_index_next = 0;
                bin_axial_next   = bin_axial; 
                bin_azimuthal_next    = bin_azimuthal; 
                if ( bin_azimuthal == (N_azimuthal-1) ) {
                     bin_axial_next += 1;
                     bin_azimuthal_next   = 0; 
                }
                else
                     bin_azimuthal_next   += 1;
                offsets_index_next = bin_azimuthal_next*N_axial + bin_axial_next;
                N_locations_plane = offsets_matrix[offsets_index_next] - offset; 
            } 

            //fprintf(stdout,"                 N_locations_plane: %d \n",N_locations_plane ); 
            
            // 2) extract location and counts for each active location 
            for (i=0; i<N_locations_plane; i++) { 
                bin_u = locations[(offset+i)*3+0]; 
                bin_v = locations[(offset+i)*3+1]; 
                bin_w = locations[(offset+i)*3+2]; // FIXME: unutilised, for TOF
//fprintf(stdout,"N_locations_plane: %d    offset: %d    bin_u: %d    bin_v: %d    bin_w: %d \n",N_locations_plane,offset,bin_u, bin_v, bin_w );                
                // make sure that bin_u and bin_v are not larger than N_u and N_v: 
                if (bin_u >= N_u || bin_v >= N_v) { 
                    fprintf(stderr, "Unhandled error while decompressing: decompressed detector bin index points out of the detector plane. \n"); 
                    fprintf(stderr, "bin_u: %d   N_u: %d   bin_v: %d   N_v: %d \n",bin_u,N_u,bin_v,N_v);  
//                    return STATUS_UNHANDLED_ERROR; 
                } 
                else {
                    projection[bin_axial*(N_u*N_v*N_azimuthal) + bin_azimuthal*(N_u*N_v) +  bin_u*(N_v) + bin_v] = counts[offset+i]; 
                    //projection[bin_azimuthal*(N_u*N_v*N_axial) + bin_axial*(N_u*N_v) + bin_u*(N_v) + bin_v] = counts[offset+i]; 
                }
            }
        }
    }

    return status; 
}


    
//Compress a static projection
unsigned int pet_compress_projection(unsigned int *n_locations, unsigned int *n_axial, unsigned int *n_azimuthal, unsigned int *n_u, unsigned int *n_v, int *offsets_matrix, float *counts, unsigned short *locations, float *projection) 
{
    int status=STATUS_SUCCESS; 

    unsigned int N_locations = *n_locations; 
    unsigned int N_axial     = *n_axial; 
    unsigned int N_azimuthal = * n_azimuthal; 
    unsigned int N_u = *n_u; 
    unsigned int N_v = *n_v; 
    
    unsigned int bin_axial; 
    unsigned int bin_azimuthal; 
    unsigned int bin_u, bin_v, bin_w, i; 
    unsigned int offset; 
    unsigned int N_locations_plane; 

    unsigned int offsets_index; 
    unsigned int offsets_index_next; 
    unsigned int bin_azimuthal_next; 
    unsigned int bin_axial_next; 
    unsigned int index; 
    unsigned int N1 = N_u*N_v*N_azimuthal; //precompute for speed
    unsigned int N2 = N_u*N_v;      //precompute for speed
    
    // for each detector plane 
    for (bin_axial=0; bin_axial < N_axial; bin_axial++) { 
        for (bin_azimuthal=0; bin_azimuthal < N_azimuthal; bin_azimuthal++) { 
            offsets_index = bin_azimuthal*N_axial + bin_axial;
            offset = offsets_matrix[offsets_index]; 
            // the following compresses the data for a single detector plane: 
            // 1) determine the offset of the data and the number of active locations 
            //fprintf(stdout,"Uncompressing:  axial: %d   azimuthal: %d    offset: %d \n",bin_axial, bin_azimuthal, offset); 
            if ((bin_azimuthal==(N_azimuthal-1)) && (bin_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                offsets_index_next = 0;
                bin_axial_next   = bin_axial; 
                bin_azimuthal_next    = bin_azimuthal; 
                if ( bin_azimuthal == (N_azimuthal-1) ) {
                     bin_axial_next += 1;
                     bin_azimuthal_next   = 0; 
                }
                else
                     bin_azimuthal_next   += 1;
                offsets_index_next = bin_azimuthal_next*N_axial + bin_axial_next;
                N_locations_plane  = offsets_matrix[offsets_index_next] - offset; 
            } 
            
            // 2) extract location and counts for each active location 
            for (i=0; i<N_locations_plane; i++) { 
                bin_u = locations[(offset+i)*3+0]; 
                bin_v = locations[(offset+i)*3+1]; 
                bin_w = locations[(offset+i)*3+2]; // FIXME: unutilised, for TOF or energy
               
                index = bin_axial*N1 + bin_azimuthal*N2 +  bin_u*(N_v) + bin_v; 
                // make sure that bin_u and bin_v are not larger than N_u and N_v: 
                if (bin_u >= N_u || bin_v >= N_v) { 
                    fprintf(stderr, "Unhandled error while compressing: decompressed detector bin index points out of the detector plane. \n"); 
                    fprintf(stderr, "bin_u: %d   N_u: %d   bin_v: %d   N_v: %d \n",bin_u,N_u,bin_v,N_v); 
                    return STATUS_UNHANDLED_ERROR; 
                } 
                else {
                    if (offset+i >= N_locations) {
                        fprintf(stderr, "Unhandled error while compressing: decompressed detector bin index points out of the detector plane. \n"); 
                        fprintf(stderr, "offset: %d   N_locations: %d \n",offset+i,N_locations);
                        return STATUS_UNHANDLED_ERROR; 
                    }
                    else
                        projection[offset+i] = counts[index]; 
                }
            }
        }
    }

    return status; 
}



//Create offsets matrix and locations array for full sampling. 
unsigned int pet_initialize_compression_structure(unsigned int *N_axial, unsigned int *N_azimuthal, unsigned int *N_u, unsigned int *N_v, int *offsets_matrix, unsigned short *locations )
{
    int status = STATUS_SUCCESS; 
    
    unsigned int n_axial         = *N_axial; 
    unsigned int n_azimuthal     = *N_azimuthal; 
    unsigned int n_u             = *N_u; 
    unsigned int n_v             = *N_v; 
    unsigned int bin_axial, bin_azimuthal, u, v;
    unsigned int offset = 0; 

    unsigned int offsets_index; 

    // for each detector plane 
    for (bin_axial=0; bin_axial < n_axial; bin_axial++) { 
        for (bin_azimuthal=0; bin_azimuthal < n_azimuthal; bin_azimuthal++) { 
            offsets_index = bin_azimuthal*n_axial + bin_axial;
            offsets_matrix[offsets_index] = offset; 
            // fill up locations array for detector plane [bin_azimuthal,bin_axial]: 
            int i = 0;
            for (u=0; u < n_u; u++) {
                for (v=0; v < n_v; v++) {
                    locations[(offset+i)*3+0] = u;
                    locations[(offset+i)*3+1] = v;
                    locations[(offset+i)*3+2] = 0; 
                    i += 1;
                }
            }
            // update offset: 
            offset += n_u*n_v; 
        }
    }
 
    return status; 
}


unsigned int pet_compress_projection_array(unsigned int *_N_ax, unsigned int *_N_az, unsigned int *_N_u, unsigned int *_N_v, float* _threshold, unsigned int *_N_active, bool *active, float *data, int *offsets, unsigned short *locations, float *counts)
{
    int status = STATUS_SUCCESS; 
    
    unsigned int N_active_locations = _N_active[0]; 
    unsigned int N_ax = _N_ax[0];
    unsigned int N_az = _N_az[0];
    unsigned int N_u = _N_u[0];
    unsigned int N_v = _N_v[0];
    float threshold = _threshold[0];
    unsigned int counter = 0; 
    unsigned int bin_ax, bin_az, bin_u, bin_v, index; 
    unsigned int N1 = N_u*N_v*N_az; //precompute for speed
    unsigned int N2 = N_u*N_v;      //precompute for speed

    //fprintf(stdout,"Events per time bin: %d \n",N_events_time_bin); 
    for (bin_ax=0; bin_ax<N_ax; bin_ax++) {
        for (bin_az=0; bin_az<N_az; bin_az++) { 
            // new plane: update offsets matrix: 
            offsets[bin_az*N_ax + bin_ax] = counter;  
            for (bin_u=0; bin_u<N_u; bin_u++) {
                for (bin_v=0; bin_v<N_v; bin_v++) {
                    index = bin_ax*N1 + bin_az*N2 +  bin_u*(N_v) + bin_v; 
                    if (active[index]) {
                    //if (data[index] > threshold) {
                        // append event: add to the array of locations 
                        if (counter < N_active_locations) {
                            counts[counter] = data[index]; 
                            locations[counter*3+0] = bin_u; 
                            locations[counter*3+1] = bin_v; 
                            locations[counter*3+2] = 0; 
                        }
                        else {
                            fprintf(stderr,"Unexpected error: counter >= N_active_locations. \n");
                            return STATUS_UNHANDLED_ERROR; 
                        }
                        counter += 1; 
                    }
                }
            }
        }
    }
    //fprintf(stdout,"Compression: counter: %d, N_locations: %d \n", counter, N_active_locations);
    //fprintf(stdout,"Compression: N_ax: %d, N_az: %d, N_u: %d, N_v: %d  \n", N_ax, N_az, N_u, N_v);
    return status;
}



unsigned int pet_get_subset_sparsity(unsigned int *_N_ax, unsigned int *_N_az, unsigned int *_N_u, unsigned int *_N_v, unsigned int *_N_locations, int *offsets, unsigned short *locations, unsigned int *subsets_matrix, int *offsets_sub, unsigned short *locations_sub)
{
    int status=STATUS_SUCCESS; 

    unsigned int N_locations = *_N_locations; 
    unsigned int N_axial     = *_N_ax; 
    unsigned int N_azimuthal = *_N_az; 
    unsigned int N_u = *_N_u; 
    unsigned int N_v = *_N_v; 
    
    unsigned int bin_axial; 
    unsigned int bin_azimuthal; 
    unsigned int bin_u, bin_v, bin_w, jj; 
    unsigned int offset; 
    unsigned int N_locations_plane; 

    unsigned int offsets_index; 
    unsigned int offsets_index_next; 
    unsigned int bin_azimuthal_next; 
    unsigned int bin_axial_next; 

    unsigned int N1 = N_u*N_v*N_azimuthal; //precompute for speed
    unsigned int N2 = N_u*N_v;      //precompute for speed
    unsigned int ii = 0; 
    unsigned int N_subsets = 0; 
    for (bin_axial=0; bin_axial < N_axial; bin_axial++) { 
        for (bin_azimuthal=0; bin_azimuthal < N_azimuthal; bin_azimuthal++) {
            offsets_index = bin_azimuthal*N_axial + bin_axial;
            if (subsets_matrix[offsets_index]) {  
                N_subsets ++; 
            }
        }
    }
    
    // for each detector plane 
    for (bin_axial=0; bin_axial < N_axial; bin_axial++) { 
        for (bin_azimuthal=0; bin_azimuthal < N_azimuthal; bin_azimuthal++) { 
            offsets_index = bin_azimuthal*N_axial + bin_axial;
            offset = offsets[offsets_index]; 
            // the following compresses the data for a single detector plane: 
            // 1) determine the offset of the data and the number of active locations 
            //fprintf(stdout,"Uncompressing:  axial: %d   azimuthal: %d    offset: %d \n",bin_axial, bin_azimuthal, offset); 
            if ((bin_azimuthal==(N_azimuthal-1)) && (bin_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                offsets_index_next = 0;
                bin_axial_next   = bin_axial; 
                bin_azimuthal_next    = bin_azimuthal; 
                if ( bin_azimuthal == (N_azimuthal-1) ) {
                     bin_axial_next += 1;
                     bin_azimuthal_next   = 0; 
                }
                else
                     bin_azimuthal_next   += 1;
                offsets_index_next = bin_azimuthal_next*N_axial + bin_axial_next;
                N_locations_plane  = offsets[offsets_index_next] - offset; 
            } 

            
            if (subsets_matrix[offsets_index]) { 
                ii++;
                if (ii<N_subsets) {
                    offsets_sub[ii] = offsets_sub[ii-1]+N_locations_plane; 
                }
                // FIXME: use memcopy
                for(jj=0; jj<N_locations_plane; jj++) {
                    locations_sub[(offsets_sub[ii-1]+jj)*3+0] = locations[(offset+jj)*3+0]; 
                    locations_sub[(offsets_sub[ii-1]+jj)*3+1] = locations[(offset+jj)*3+1]; 
                    locations_sub[(offsets_sub[ii-1]+jj)*3+2] = locations[(offset+jj)*3+2]; 
                }
            }            
        }
    }
    return status; 
}



unsigned int pet_get_subset_projection_array(unsigned int *_N_ax, unsigned int *_N_az, unsigned int *_N_u, unsigned int *_N_v, unsigned int *_N_locations, float *data, int *offsets, unsigned short *locations, int *offsets_sub, unsigned short *locations_sub, unsigned int *subsets_matrix, float *data_sub)
{
    int status=STATUS_SUCCESS; 

    unsigned int N_locations = *_N_locations; 
    unsigned int N_axial     = *_N_ax; 
    unsigned int N_azimuthal = *_N_az; 
    unsigned int N_u = *_N_u; 
    unsigned int N_v = *_N_v; 
    
    unsigned int bin_axial; 
    unsigned int bin_azimuthal; 
    unsigned int bin_u, bin_v, bin_w, jj; 
    unsigned int offset; 
    unsigned int N_locations_plane; 

    unsigned int offsets_index; 
    unsigned int offsets_index_next; 
    unsigned int bin_azimuthal_next; 
    unsigned int bin_axial_next; 
    unsigned int index; 

    unsigned int N1 = N_u*N_v*N_azimuthal; //precompute for speed
    unsigned int N2 = N_u*N_v;      //precompute for speed
    unsigned int ii = 0; 
    unsigned int N_subsets = 0; 
    for (bin_axial=0; bin_axial < N_axial; bin_axial++) { 
        for (bin_azimuthal=0; bin_azimuthal < N_azimuthal; bin_azimuthal++) {
            offsets_index = bin_azimuthal*N_axial + bin_axial;
            if (subsets_matrix[offsets_index]) {  
                N_subsets ++; 
            }
        }
    }
    
    // for each detector plane 
    for (bin_azimuthal=0; bin_azimuthal < N_azimuthal; bin_azimuthal++) { 
        for (bin_axial=0; bin_axial < N_axial; bin_axial++) { 
            offsets_index = bin_azimuthal*N_axial + bin_axial;
            offset = offsets[offsets_index]; 
            // the following compresses the data for a single detector plane: 
            // 1) determine the offset of the data and the number of active locations 
            //fprintf(stdout,"Uncompressing:  axial: %d   azimuthal: %d    offset: %d \n",bin_axial, bin_azimuthal, offset); 
            if ((bin_azimuthal==(N_azimuthal-1)) && (bin_axial==(N_axial-1))) { 
                N_locations_plane = N_locations - offset; 
            }
            else {
                offsets_index_next = 0;
                bin_axial_next   = bin_axial; 
                bin_azimuthal_next    = bin_azimuthal; 
                if ( bin_azimuthal == (N_azimuthal-1) ) {
                     bin_axial_next += 1;
                     bin_azimuthal_next   = 0; 
                }
                else
                     bin_azimuthal_next   += 1;
                offsets_index_next = bin_azimuthal_next*N_axial + bin_axial_next;
                N_locations_plane  = offsets[offsets_index_next] - offset; 
            } 

            
            if (subsets_matrix[offsets_index]) { 
                ii++;
                if (ii<N_subsets) {
                    offsets_sub[ii] = offsets_sub[ii-1]+N_locations_plane; 
                }
                // FIXME: use memcopy
                for(jj=0; jj<N_locations_plane; jj++) {
                    locations_sub[(offsets_sub[ii-1]+jj)*3+0] = locations[(offset+jj)*3+0]; 
                    locations_sub[(offsets_sub[ii-1]+jj)*3+1] = locations[(offset+jj)*3+1]; 
                    locations_sub[(offsets_sub[ii-1]+jj)*3+2] = locations[(offset+jj)*3+2]; 
                    
                    index = bin_axial*N1 + bin_azimuthal*N2 +  locations[(offset+jj)*3+0]*(N_v) + locations[(offset+jj)*3+1]; 
                    
                    data_sub[offsets_sub[ii-1]+jj] = data[index]; 
                    //data_sub[offsets_sub[ii-1]+jj] = data[offset+jj]; 
                }
            }            
        }
    }
    return status; 
}

