/*
 *  _tt.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_tt.h"
#include "_et_common.h"

////////////////////////////////////////////////////////////////////////////////
// tt_project_ray
////////////////////////////////////////////////////////////////////////////////
int tt_project_ray(VolumeType h_volume[], u_int_3 volume_voxels, float out_projections[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_position[], float_3 volume_size, float t_step, int use_gpu)
{
#ifdef _USE_CUDA
    dim3 blockSize(16, 16);
    bool linearFiltering = true;
    dim3 gridSize;

    float  invViewMatrix[12];
    float3 sourcePosition;
    float3 volumeSize;
    float *d_output;
    float *d_output_proj;

    //cuInit(0);

    cudaExtent vsize = make_cudaExtent(volume_voxels.x,volume_voxels.y,volume_voxels.z);
//fprintf(stderr, "\n %d, %d, %d",volume_voxels.x,volume_voxels.y,volume_voxels.z);
    initCuda(h_volume, vsize);
    setTextureFilterMode(linearFiltering);
    gridSize = dim3(iDivUp(detector_pixels.w, blockSize.x), iDivUp(detector_pixels.h, blockSize.y));

    //Allocate memory for projections on the device
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_output, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float) ));
    CUDA_SAFE_CALL(cudaMemset((void *)d_output,0, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float) ));

 //   struct timeval start_time; gettimeofday( &start_time, 0);
 //   struct timeval t_time;
 //   float elapsed_time;

    volumeSize.x = volume_size.x;
    volumeSize.y = volume_size.y;
    volumeSize.z = volume_size.z;

    for (int proj=0;proj<n_projections;proj++)
    {
        //define invViewMatrix (position of detector) and position of source
        set_inViewMatrix(invViewMatrix, detector_scale[proj], detector_transl[proj], detector_rotat[proj]);

        sourcePosition.x = source_position[proj].x;
        sourcePosition.y = source_position[proj].y;
        sourcePosition.z = source_position[proj].z;

        //project
        copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);
        d_output_proj = (float*) d_output + proj * detector_pixels.w * detector_pixels.h;
        tt_line_project_ray_gpu(gridSize, blockSize, d_output_proj, sourcePosition, volumeSize, detector_pixels.w, detector_pixels.h, t_step);
    }

 //   gettimeofday( &t_time, 0);
 //   elapsed_time = (float) (1000.0 * ( t_time.tv_sec - start_time.tv_sec) + (0.001 * (t_time.tv_usec - start_time.tv_usec)) );
 //   fprintf(stderr,"\nTime per projection %d %d %d -> %d %d: %f ms",volume_voxels.x,volume_voxels.y,volume_voxels.z,detector_pixels.w,detector_pixels.h,elapsed_time/n_projections);

    //Copy result back to host
    CUDA_SAFE_CALL(cudaMemcpy(out_projections, d_output, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float), cudaMemcpyDeviceToHost));

    //Clean up
    CUDA_SAFE_CALL(cudaFree(d_output));
    freeCudaBuffers();
    cudaThreadExit();
    return 0;
#else
    return 1;
#endif
}



////////////////////////////////////////////////////////////////////////////////
// tt_backproject_ray
////////////////////////////////////////////////////////////////////////////////

int tt_backproject_ray(float h_projections[], u_int_2 detector_pixels, u_int n_projections, float out_backprojection[], float_2 detector_size[], float_3 detector_transl[], float_3 detector_rotat[], float_3 source_pos[], u_int_3 volume_voxels, float_3 volume_size, float t_step, int interpolation, int use_gpu)
{
#ifdef _USE_CUDA

    float  invViewMatrix[12];
    float3 sourcePosition;
    float3 volumeSize;
    uint3  volumeVoxels;
    uint2  detectorPixels;

    float *d_projections;
    float *d_output;
    float *d_current_projection;

    dim3 blockSize(16, 16);
    dim3 gridSize;
    if (use_gpu)
    {
      //cuInit(0);
      CUDA_SAFE_CALL(cudaMalloc((void **)&d_projections, detector_pixels.w*detector_pixels.h*n_projections*sizeof(float) ));
      CUDA_SAFE_CALL(cudaMemcpy(d_projections, h_projections, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float), cudaMemcpyHostToDevice));
      gridSize = dim3(iDivUp(detector_pixels.w, blockSize.x), iDivUp(detector_pixels.h, blockSize.y));

      //Allocate memory for backprojection on the device
      CUDA_SAFE_CALL(cudaMalloc((void **)&d_output, volume_voxels.x*volume_voxels.y*volume_voxels.z*sizeof(float) ));
      CUDA_SAFE_CALL(cudaMemset((void *)d_output,0, volume_voxels.x*volume_voxels.y*volume_voxels.z*sizeof(float) ));
    }

//    struct timeval start_time; gettimeofday( &start_time, 0);
//    struct timeval t_time;
//    float elapsed_time;

    volumeSize.x = volume_size.x;
    volumeSize.y = volume_size.y;
    volumeSize.z = volume_size.z;
    volumeVoxels.x = volume_voxels.x;
    volumeVoxels.y = volume_voxels.y;
    volumeVoxels.z = volume_voxels.z;
    detectorPixels.x = detector_pixels.w;
    detectorPixels.y = detector_pixels.h;

//for (int proj=49;proj<50;proj++) 
    for (int proj=0;proj<n_projections;proj++)
    {
        //define invViewMatrix (position of detector) and position of source
//        fprintf(stderr,"\n----> %d / %d - %f %f - %f %f %f - %f %f %f ",proj,n_projections,detector_size[proj].x,detector_size[proj].y,detector_transl[proj].x,detector_transl[proj].y,detector_transl[proj].z,detector_rotat[proj].x,detector_rotat[proj].y,detector_rotat[proj].z);
        set_inViewMatrix(invViewMatrix, detector_size[proj], detector_transl[proj], detector_rotat[proj]);
        sourcePosition.x = source_pos[proj].x;
        sourcePosition.y = source_pos[proj].y;
        sourcePosition.z = source_pos[proj].z;
        //backproject

        if (use_gpu)
        {
            copyInvViewMatrix_bk(invViewMatrix, sizeof(float4)*3);
            d_current_projection = (float*) d_projections + proj * detector_pixels.w * detector_pixels.h;
            tt_line_backproject_ray_gpu(gridSize, blockSize, d_current_projection, d_output, detectorPixels, sourcePosition, volumeVoxels, volumeSize, t_step, interpolation);
 //           fprintf(stderr,"\n => %d %d %d - %d %d %d - %d %d - %f %f %f - %d %d %d - %f %f %f - %f",gridSize.x,gridSize.y,gridSize.z, blockSize.x,blockSize.y,blockSize.z, detectorPixels.x,detectorPixels.y, sourcePosition.x,sourcePosition.y,sourcePosition.z, volumeVoxels.x,volumeVoxels.y,volumeVoxels.z, volumeSize.x,volumeSize.y,volumeSize.z, t_step);
        }
        else
        {

            float *current_projection = (float*) h_projections + proj * detector_pixels.w * detector_pixels.h;
//            fprintf(stderr,"\nBack-projection %d/%d",proj+1,n_projections);
//            tt_line_backproject_ray_cpu(out_backprojection, current_projection, invViewMatrix, detectorPixels, sourcePosition, volumeVoxels, volumeSize, t_step, interpolation);

        }

    }
//    gettimeofday( &t_time, 0);
//    elapsed_time = (float) (1000.0 * ( t_time.tv_sec - start_time.tv_sec) + (0.001 * (t_time.tv_usec - start_time.tv_usec)) );
//    fprintf(stderr,"\nTime per back-projection %d x [%d %d] -> [%d %d %d]: %f ms",n_projections,detector_pixels.w,detector_pixels.h,volume_voxels.x,volume_voxels.y,volume_voxels.z,elapsed_time/n_projections);

//cudaMemset(d_output, 100, volume_voxels.x*volume_voxels.y*volume_voxels.z*sizeof(float) );

    //Copy result back to host

    if(use_gpu)
    {
        CUDA_SAFE_CALL(cudaMemcpy(out_backprojection, d_output, volume_voxels.x*volume_voxels.y*volume_voxels.z*sizeof(float), cudaMemcpyDeviceToHost));

        //Clean up
        CUDA_SAFE_CALL(cudaFree(d_output));
        CUDA_SAFE_CALL(cudaFree(d_projections));

        cudaThreadExit();
    }
    return 0;
#else
    return 1;
#endif
}






////////////////////////////////////////////////////////////////////////////////
// tt_project
////////////////////////////////////////////////////////////////////////////////

int tt_project_perspective(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background)
{
    int status = 1;
    status = 0;
    return status;
}


////////////////////////////////////////////////////////////////////////////////
// tt_project_gpu
////////////////////////////////////////////////////////////////////////////////

#ifdef _USE_CUDA
int tt_project_perspective_gpu(nifti_image *attenuation, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background)
{
/*	// initialise the cuda arrays //
        cudaArray *attenuationArray_d;            //stores input attenuation coefficients, makes use of fetch unit
	float     *projectionArray_d;             //stores projections (output)
        float     *perspectiveAttenuationArray_d; //stores attenuation coefficients aligned with current camera
	float4    *positionFieldImageArray_d;     //stores position field for rotation of attenuation
	int       *mask_d;                        //binary mask that defines active voxels (typically all active)
        float     *psfArray_d;                    //stores point spread function
        int       psf_size[3];
        int       image_size[3];

	// Allocate arrays on the device //
	if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuation->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&projectionArray_d, projectionImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&perspectiveAttenuationArray_d, attenuation->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, attenuation->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, attenuation->dim)) return 1;

	// Transfer data from the host to the device //
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuation)) return 1;
	int *mask_h=(int *)malloc(attenuation->nvox*sizeof(int));
	for(int i=0; i<attenuation->nvox; i++) mask_h[i]=i;
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, attenuation->nvox*sizeof(int), cudaMemcpyHostToDevice));
	free(mask_h);

	// Alloc transformation matrix //
	mat44 *transformationMatrix = (mat44 *)calloc(1,sizeof(mat44));

        // Allocate and initialize kernel for DDPSF //
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            image_size[0] = attenuation->dim[1];
            image_size[1] = attenuation->dim[2];
            image_size[2] = attenuation->dim[3];
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            }

	for(unsigned int projection=0; projection<n_projections; projection++){

		// Compute deformation field //

//		et_create_rotation_matrix(transformationMatrix, 0,0,0,0,0,0);

//		reg_affine_positionField_gpu(	transformationMatrix,
//						attenuation,
//						&positionFieldImageArray_d);

                tt_perspective_positionField_gpu( attenuation,
                                                  image_origin, 
                                                  detector_origin,
                                                  detector_shape,
                                                  &positionFieldImageArray_d);
fprintf(stderr, "iO: %f %f %f\n dO: %f %f %f\n dS: %f %f \n",image_origin[0],image_origin[1],image_origin[2],detector_origin[0],detector_origin[1],detector_origin[2], detector_shape[0], detector_shape[1]);

		// Resample attenuation image //
		reg_resampleSourceImage_gpu(	attenuation,
						attenuation,
						&perspectiveAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuation->nvox,
						background);

		// Correct for non-uniform voxel size //

		// Integrate along lines //
                et_line_integral_gpu(		&perspectiveAttenuationArray_d,
						&projectionArray_d,
						projection,
						attenuation);
	}

	// Transfer result back to host //
	if(cudaCommon_transferFromDeviceToNifti(projectionImage, &projectionArray_d)) return 1;

	//Free//
	cudaCommon_free(&attenuationArray_d);
	cudaCommon_free((void **)&perspectiveAttenuationArray_d);
	cudaCommon_free((void **)&projectionArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
	free(transformationMatrix);
        if (psfImage != NULL)
            cudaCommon_free((void **)&psfArray_d);
	return 0; 
*/
return 1;
}
#endif


