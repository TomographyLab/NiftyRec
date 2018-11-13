
#include "_reg.h"

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_bspline.h"
#include "_reg_mutualinformation.h"
#include "_reg_ssd.h"
#include "_reg_tools.h"
#include "float.h"

#ifdef _USE_CUDA
	#include "_reg_cudaCommon.h"
	#include "_reg_resampling_gpu.h"
	#include "_reg_affineTransformation_gpu.h"
	#include "_reg_bspline_gpu.h"
	#include "_reg_mutualinformation_gpu.h"
	#include "_reg_tools_gpu.h"
#endif

#define JH_PW_APPROX 2


int reg_gradient_NMI_nodes(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointImage, nifti_image* nodeGradientImage, int binning)
{
	int status = 1;
	return status;
}

#ifdef _USE_CUDA
int reg_gradient_NMI_nodes_gpu(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointImage, nifti_image* nodeGradientImage, int binning)
{
	float *targetImageArray_d=NULL;
	cudaArray *sourceImageArray_d=NULL;
	float4 *controlPointImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	int *targetMask_d=NULL;
	float4 *resultGradientArray_d=NULL;
	float4 *voxelNMIGradientArray_d=NULL;
	float4 *nodeNMIGradientArray_d=NULL;
	float *resultImageArray_d=NULL;
	float *logJointHistogram_d=NULL;

        double *probaJointHistogram = (double *)malloc(binning*(binning+2)*sizeof(double));
        double *logJointHistogram = (double *)malloc(binning*(binning+2)*sizeof(double));
        double *entropies = (double *)malloc(4*sizeof(double));
//        float *probaJointHistogram = (float *)malloc(binning*(binning+2)*sizeof(float));
//        float *logJointHistogram = (float *)malloc(binning*(binning+2)*sizeof(float));
//	float *entropies = (float *)malloc(4*sizeof(float));
//	double *entropies_double = (double *)malloc(4*sizeof(double));	

	int smoothingRadius[3];
	smoothingRadius[0] = (int)floor( 2.0*controlPointImage->dx/targetImage->dx );
	smoothingRadius[1] = (int)floor( 2.0*controlPointImage->dy/targetImage->dy );
	smoothingRadius[2] = (int)floor( 2.0*controlPointImage->dz/targetImage->dz );

	int activeVoxelNumber = targetImage->nvox;

	if(cudaCommon_allocateArrayToDevice<float>(&targetImageArray_d, targetImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&targetImageArray_d, targetImage)) return 1;
fprintf(stderr,"\n1.\n");
	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;
fprintf(stderr,"\n1..\n");
	if(cudaCommon_allocateArrayToDevice<float>(&resultImageArray_d, targetImage->dim)) return 1;
fprintf(stderr,"\n1...\n");
        fprintf(stderr,"\nControl point image dim: %d %d %d %d %d %d\n",controlPointImage->dim[0],controlPointImage->dim[1],controlPointImage->dim[2],controlPointImage->dim[3],controlPointImage->dim[4],controlPointImage->dim[5]);
	if(cudaCommon_allocateArrayToDevice<float4>(&controlPointImageArray_d, controlPointImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float4>(&controlPointImageArray_d,controlPointImage)) return 1;
fprintf(stderr,"\n1....\n");
	CUDA_SAFE_CALL(cudaMalloc((void **)&logJointHistogram_d, binning*(binning+2)*sizeof(float)));
	if(cudaCommon_allocateArrayToDevice<float4>(&voxelNMIGradientArray_d, sourceImage->dim)) return 1; //result
	if(cudaCommon_allocateArrayToDevice<float4>(&nodeNMIGradientArray_d, controlPointImage->dim)) return 1;
fprintf(stderr,"\n1.....\n");
CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
int *mask_h=(int *)malloc(activeVoxelNumber*sizeof(int));
for(int i=0; i<activeVoxelNumber; i++) mask_h[i]=i;
CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, mask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));
fprintf(stderr,"\n1......\n");
	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&resultGradientArray_d, activeVoxelNumber*sizeof(float4)));
fprintf(stderr,"\n3\n");
	nifti_image *tempImage=NULL;
	tempImage = nifti_copy_nim_info(sourceImage);
	tempImage->datatype = NIFTI_TYPE_FLOAT32;
	tempImage->nbyper = sizeof(float);
	tempImage->data = (float *)malloc(sourceImage->nvox*sizeof(float));
fprintf(stderr,"\n4\n");
	// generate the position field //
	reg_bspline_gpu(			controlPointImage,
						targetImage,
						&controlPointImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber);
fprintf(stderr,"\n5\n");
	// Resample the source image //
	reg_resampleSourceImage_gpu(		targetImage, //result
						sourceImage,
						&resultImageArray_d,
						&sourceImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber,
						0,
						0);
fprintf(stderr,"\n6\n");
	// The result image is transfered back to the host //
	if(cudaCommon_transferFromDeviceToNifti(tempImage, &resultImageArray_d)) return 1;

//nifti_set_filenames(targetImage, "/home/spedemon/Desktop/test/test_target.nii", 0, 0);
//nifti_image_write(targetImage);
//nifti_set_filenames(sourceImage, "/home/spedemon/Desktop/test/test_source.nii", 0, 0);
//nifti_image_write(sourceImage);
//nifti_set_filenames(tempImage, "/home/spedemon/Desktop/test/test_temp.nii", 0, 0);
//nifti_image_write(tempImage);
//fprintf(stderr,"Binning: %d \n",binning);

//for(int i=0; i<binning*(binning+2); i++)
//    fprintf(stderr,"%f ",logJointHistogram[i]);
//fprintf(stderr,"\n");
//for(int i=0; i<binning*(binning+2); i++)
//    fprintf(stderr,"%f ",probaJointHistogram[i]);
//fprintf(stderr,"\n");
fprintf(stderr,"Types: %d %d %d \n",targetImage->datatype, sourceImage->datatype, NIFTI_TYPE_FLOAT32);
//fprintf(stderr,"targetImage: %d %d %d %d \n",targetImage->nx,targetImage->ny,targetImage->nz,targetImage->nvox);
//fprintf(stderr,"tempImage:   %d %d %d %d \n",tempImage->nx,tempImage->ny,tempImage->nz,tempImage->nvox);


	/* Create joint histogram and transfer to the device*/
	reg_getEntropies<double>(		targetImage,
						tempImage, //resultImage
						2,
						binning,
						probaJointHistogram,
						logJointHistogram,
						entropies,
						mask_h);

//for(int i=0; i<binning*(binning+2); i++)
//    fprintf(stderr,"%f ",logJointHistogram[i]);
//fprintf(stderr,"\n");
//for(int i=0; i<binning*(binning+2); i++)
//    fprintf(stderr,"%f ",probaJointHistogram[i]);
//fprintf(stderr,"\n");

//	entropies_double[0]=entropies[0];entropies_double[1]=entropies[1];entropies_double[2]=entropies[2];entropies_double[3]=entropies[3];
fprintf(stderr,"\nBinning: %d  Entropies: %f %f %f %f",binning, entropies[0],entropies[1],entropies[2],entropies[3]);


	// Tranfer the histogram to the device //
	CUDA_SAFE_CALL(cudaMemcpy(logJointHistogram_d, logJointHistogram, binning*(binning+2)*sizeof(float), cudaMemcpyHostToDevice));

	// NMI Gradient //
	reg_getSourceImageGradient_gpu(						targetImage,
										sourceImage,
										&sourceImageArray_d,
										&positionFieldImageArray_d,
										&resultGradientArray_d,
										activeVoxelNumber,
										0);

	reg_getVoxelBasedNMIGradientUsingPW_gpu(				targetImage,
										tempImage, 
										&targetImageArray_d,     // V
										&resultImageArray_d,     // V
										&resultGradientArray_d,  // V
										&logJointHistogram_d,    // V
										&voxelNMIGradientArray_d,
										&targetMask_d,           // V
										activeVoxelNumber,       //
										entropies,        //
										binning);                //

fprintf(stderr,"\ngradientImage: %d %d %d\n",nodeGradientImage->nx,nodeGradientImage->ny,nodeGradientImage->nz);
fprintf(stderr,"sourceImage: %d %d %d\n",sourceImage->nx,sourceImage->ny,sourceImage->nz);
fprintf(stderr,"Smoothing radius: %d %d %d\n", smoothingRadius[0],smoothingRadius[1],smoothingRadius[2]);
fprintf(stderr,"Image size: %d %d %d\n", sourceImage->nx, sourceImage->ny, sourceImage->nz);

	reg_smoothImageForCubicSpline_gpu(					tempImage,
										&voxelNMIGradientArray_d,
										smoothingRadius);

//	reg_voxelCentric2NodeCentric_gpu(					tempImage,
//										controlPointImage,
//										&voxelNMIGradientArray_d,
//										&nodeNMIGradientArray_d);

        reg_voxelCentric2NodeCentric_gpu(                                       targetImage,
                                                                                controlPointImage,
                                                                                &voxelNMIGradientArray_d,
                                                                                &nodeNMIGradientArray_d);

//				reg_convertNMIGradientFromVoxelToRealSpace_gpu( sourceMatrix_xyz,
//										                        controlPointImage,
//										                        &nodeNMIGradientArray_d);


fprintf(stderr,"Voxel number: %d\n",activeVoxelNumber);
//fprintf(stderr,"Nodes size:   %d\n",nodeGradientImage->nvox);

//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, nodeNMIGradientArray_d, nodeGradientImage->nvox*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, logJointHistogram_d, binning*(binning+2)*sizeof(float), cudaMemcpyDeviceToHost));
CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, voxelNMIGradientArray_d, activeVoxelNumber*4*sizeof(float), cudaMemcpyDeviceToHost));

//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, resultGradientArray_d, activeVoxelNumber*4*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, positionFieldImageArray_d, activeVoxelNumber*4*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, targetMask_d, activeVoxelNumber*1*sizeof(float), cudaMemcpyDeviceToHost));

	cudaCommon_free( (void **)&targetImageArray_d );
	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&controlPointImageArray_d );
	cudaCommon_free( (void **)&resultImageArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	CUDA_SAFE_CALL(cudaFree(targetMask_d));
	cudaCommon_free((void **)&resultGradientArray_d);
	cudaCommon_free((void **)&voxelNMIGradientArray_d);
	cudaCommon_free((void **)&nodeNMIGradientArray_d);
	cudaCommon_free((void **)&logJointHistogram_d);

        free(probaJointHistogram);
        free(logJointHistogram);
	free(entropies);
        free(mask_h);

        nifti_image_free(tempImage);

fprintf(stderr,"_reg Done \n\n");
	return 0;
}
#endif



int reg_gradient_voxel_to_nodes(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage)
{
	int status = 1;
	return status;
}


#ifdef _USE_CUDA
int reg_gradient_voxel_to_nodes_gpu(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage)
{
	int status = 1;
//	float4 *voxelNMIGradientArray_d=NULL;
//	float4 *nodeNMIGradientArray_d=NULL;

//	if(cudaCommon_allocateArrayToDevice(&voxelNMIGradientArray_d, gradientImage->dim)) return 1;
//	if(cudaCommon_allocateArrayToDevice(&nodeNMIGradientArray_d, controlPointImage->dim)) return 1;	

//	int smoothingRadius[3];
//	smoothingRadius[0] = (int)floor( 2.0*controlPointImage->dx/gradientImage->dx );
//	smoothingRadius[1] = (int)floor( 2.0*controlPointImage->dy/gradientImage->dy );
//	smoothingRadius[2] = (int)floor( 2.0*controlPointImage->dz/gradientImage->dz );

//fprintf(stderr, "Image delta:       %f %f %f \n",gradientImage->dx,gradientImage->dy,gradientImage->dz);
//fprintf(stderr, "Image size:        %d %d %d (nvox: %d) \n",gradientImage->nx,gradientImage->ny,gradientImage->nz,gradientImage->nvox);
//fprintf(stderr, "CP delta:          %f %f %f\n",controlPointImage->dx,controlPointImage->dy,controlPointImage->dz);
//fprintf(stderr, "Smoothing R:       %f \n",smoothingRadius);

//fprintf(stderr, "CP gradient delta: %f %f %f \n",cpGradientImage->dx,cpGradientImage->dy,cpGradientImage->dz);
//fprintf(stderr, "CP gradient size:  %d %d %d (nvox: %d) \n",cpGradientImage->nx,cpGradientImage->ny,cpGradientImage->nz,cpGradientImage->nvox);
//fprintf(stderr, "dim grad: %d   dim cp: %d \n",gradientImage->dim[0],controlPointImage->dim[0]);


//	if(cudaCommon_transferNiftiToArrayOnDevice<float4> (&voxelNMIGradientArray_d,gradientImage)) return 1;

//	reg_smoothImageForCubicSpline_gpu(					gradientImage,
//										&voxelNMIGradientArray_d,
//										smoothingRadius);

//nifti_set_filenames(controlPointImage, "/home/spedemon/Desktop/test/test_cp.nii", 0, 0);
//nifti_image_write(controlPointImage);

//nifti_set_filenames(gradientImage, "/home/spedemon/Desktop/test/test_gradient.nii", 0, 0);
//nifti_image_write(gradientImage);

//	reg_voxelCentric2NodeCentric_gpu(					controlPointImage,
//										gradientImage,
//										&voxelNMIGradientArray_d,
//										&nodeNMIGradientArray_d);
        reg_voxelCentric2NodeCentric(                                           cpGradientImage,
							                        gradientImage);
//	CUDA_SAFE_CALL(cudaMemcpy(cpGradientImage->data, nodeNMIGradientArray_d, cpGradientImage->nvox*sizeof(float), cudaMemcpyDeviceToHost));

//	cudaCommon_free((void **)&voxelNMIGradientArray_d);
//	cudaCommon_free((void **)&nodeNMIGradientArray_d);
	status = 0;
	return status;
}
#endif


int reg_resample_spline(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage)
{
	int status = 1;
	return status;
}


#ifdef _USE_CUDA
int reg_resample_spline_gpu(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage)
{
	int status = 1;

	cudaArray *sourceImageArray_d=NULL;
	float4 *controlPointImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	int *targetMask_d=NULL;

	float *resultImageArray_d=NULL;

	int *targetMask;
	targetMask = (int *)malloc(sourceImage->nvox*sizeof(int));
	int activeVoxelNumber=0;
	for(unsigned int i=0; i<sourceImage->nvox; i++)
		targetMask[i]=i;
	activeVoxelNumber=sourceImage->nvox;

	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;

	if(cudaCommon_allocateArrayToDevice<float>(&resultImageArray_d, sourceImage->dim)) return 1;

	if(cudaCommon_allocateArrayToDevice<float4>(&controlPointImageArray_d, controlPointImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float4>(&controlPointImageArray_d,controlPointImage)) return 1;

	// Index of the active voxel is stored
	int *targetMask_h; CUDA_SAFE_CALL(cudaMallocHost((void **)&targetMask_h, activeVoxelNumber*sizeof(int)));
	int *targetMask_h_ptr = &targetMask_h[0];
	for(unsigned int i=0;i<sourceImage->nvox;i++)
		{
		if(targetMask[i]!=-1) 
			*targetMask_h_ptr++=i;
		}
	CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, targetMask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));


	// generate the position field //
	reg_bspline_gpu(			controlPointImage,
						sourceImage,
						&controlPointImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber);

	// Resample the source image //
	reg_resampleSourceImage_gpu(		resultImage,
						sourceImage,
						&resultImageArray_d,
						&sourceImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber,
						0,
						0);

	if(cudaCommon_transferFromDeviceToNifti(resultImage, &resultImageArray_d)) return 1;

        free(targetMask);
	CUDA_SAFE_CALL(cudaFreeHost(targetMask_h));
	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&resultImageArray_d );
	cudaCommon_free( (void **)&controlPointImageArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	cudaCommon_free( (void **)&targetMask_d );

	status = 0;
	return status;
}
#endif


int reg_image_gradient(nifti_image *resultGradientImage, nifti_image *sourceImage, nifti_image *controlPointImage)
{
	int status = 1;
        float sourcePaddingValue = 0.0f;
        
        int image_size[3];
        float image_spacing[3];
        image_size[0] = sourceImage->nx; image_size[1] = sourceImage->ny; image_size[2] = sourceImage->nz;
        image_spacing[0] = sourceImage->dx; image_spacing[1] = sourceImage->dy; image_spacing[2] = sourceImage->dz;
        nifti_image *positionFieldImage = reg_initialize_deformation_field(image_size, image_spacing);
        nifti_image *resampledImage = reg_initialize_image(image_size, image_spacing);
        positionFieldImage->data = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*3*sizeof(float));
        resampledImage->data = (float*) malloc(image_size[0]*image_size[1]*image_size[2]*sizeof(float));

        /* Create mask */
	int *targetMask;
	targetMask = (int *)malloc(sourceImage->nvox*sizeof(int));
	int activeVoxelNumber=0;
	for(unsigned int i=0; i<sourceImage->nvox; i++)
		targetMask[i]=i;

	/* generate the position field */
	reg_bspline<float>(		controlPointImage,
					sourceImage, 
					positionFieldImage,
					targetMask,
					0);

	/* Resample the source image */
	reg_resampleSourceImage<float>(
					sourceImage,
					sourceImage,
					resampledImage,
					positionFieldImage,
					targetMask,
					1,
					sourcePaddingValue);

	reg_getSourceImageGradient<float>(		sourceImage,
							sourceImage,
							resultGradientImage,
							positionFieldImage,
							targetMask,
							1);
    nifti_image_free(positionFieldImage);
    nifti_image_free(resampledImage);
    free(targetMask);

    return 0;
}



#ifdef _USE_CUDA
int reg_image_gradient_gpu(nifti_image *resultGradientImage, nifti_image *sourceImage, nifti_image *controlPointImage)
{
	int status = 1;

	cudaArray *sourceImageArray_d=NULL;
	float4 *controlPointImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	float4 *resultGradientArray_d=NULL;
	int *targetMask_d=NULL;

	int *targetMask;
	targetMask = (int *)malloc(sourceImage->nvox*sizeof(int));
	int activeVoxelNumber=0;
	for(unsigned int i=0; i<sourceImage->nvox; i++)
		targetMask[i]=i;
	activeVoxelNumber=sourceImage->nvox;

	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;

	if(cudaCommon_allocateArrayToDevice<float4>(&controlPointImageArray_d, controlPointImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float4>(&controlPointImageArray_d,controlPointImage)) return 1;

	// Index of the active voxel is stored
	int *targetMask_h; CUDA_SAFE_CALL(cudaMallocHost((void **)&targetMask_h, activeVoxelNumber*sizeof(int)));
	int *targetMask_h_ptr = &targetMask_h[0];
	for(unsigned int i=0;i<sourceImage->nvox;i++)
		{
		if(targetMask[i]!=-1) 
			*targetMask_h_ptr++=i;
		}
	CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, targetMask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&resultGradientArray_d, activeVoxelNumber*sizeof(float4)));

	// generate the position field //
	reg_bspline_gpu(	controlPointImage,
						sourceImage,  
						&controlPointImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber);

        // copute gradient //
	reg_getSourceImageGradient_gpu(		sourceImage,
						sourceImage,
						&sourceImageArray_d,
						&positionFieldImageArray_d,
						&resultGradientArray_d,
						activeVoxelNumber,
						0);

	if(cudaCommon_transferFromDeviceToNifti(resultGradientImage, &resultGradientArray_d)) return 1;

	free(targetMask);
	CUDA_SAFE_CALL(cudaFreeHost(targetMask_h));
	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&controlPointImageArray_d );
	cudaCommon_free( (void **)&resultGradientArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	cudaCommon_free( (void **)&targetMask_d );

	status = 0;
	return status;
}
#endif



int reg_ssd_gradient(nifti_image *ssdGradientImage, nifti_image *targetImage, nifti_image *sourceImage, nifti_image *gradientImage, int smoothingRadius[])
{
//        float SSDValue=0.0f;
//	SSDValue = reg_getSSD<float>(			targetImage,
//							sourceImage);

/*
fprintf(stderr, "Source delta:       %f %f %f \n",sourceImage->dx,sourceImage->dy,sourceImage->dz);
fprintf(stderr, "Source size:        %d %d %d (nvox: %d) \n",sourceImage->nx,sourceImage->ny,sourceImage->nz,sourceImage->nvox);

fprintf(stderr, "Target delta:       %f %f %f \n",targetImage->dx,targetImage->dy,targetImage->dz);
fprintf(stderr, "Target size:        %d %d %d (nvox: %d) \n",targetImage->nx,targetImage->ny,targetImage->nz,targetImage->nvox);

fprintf(stderr, "Gradient delta:     %f %f %f \n",gradientImage->dx,gradientImage->dy,gradientImage->dz);
fprintf(stderr, "Gradient size:      %d %d %d (nvox: %d) \n",gradientImage->nx,gradientImage->ny,gradientImage->nz,gradientImage->nvox);

fprintf(stderr, "SSDGradient delta:  %f %f %f \n",ssdGradientImage->dx,ssdGradientImage->dy,ssdGradientImage->dz);
fprintf(stderr, "SSDGradient size:   %d %d %d (nvox: %d) \n",ssdGradientImage->nx,ssdGradientImage->ny,ssdGradientImage->nz,ssdGradientImage->nvox);
*/

	reg_getVoxelBasedSSDGradient<float>(
							targetImage,
							sourceImage,
							gradientImage,
							ssdGradientImage);
	reg_smoothImageForCubicSpline<float>(		ssdGradientImage,
							smoothingRadius);
	return 0;
}

#ifdef _USE_CUDA
int reg_ssd_gradient_gpu(nifti_image *ssdGradientImage, nifti_image *targetImage, nifti_image *sourceImage, nifti_image *gradientImage, int smoothingRadius[])
{
	int status = 1;
	return status;
}
#endif



int reg_gaussian_smooth(nifti_image *niftiImage, float sigma)
{
	bool smoothAxis[8]={true,true,true,true,true,true,true,true};
	reg_gaussianSmoothing<float>(niftiImage, sigma, smoothAxis);
	return 0;
}

#ifdef _USE_CUDA
int reg_gaussian_smooth_gpu(nifti_image *niftiImage, float sigma)
{
//	bool smoothAxis[8]={true,true,true,true,true,true,true,true};
//	reg_gaussianSmoothing<float>(niftiImage, sigma, smoothAxis);
	return 1;
}
#endif



int reg_scale_amplitude(nifti_image *niftiImage, float min_value, float max_value)
{
	reg_intensityRescale(niftiImage, min_value, max_value, -FLT_MAX, FLT_MAX);
	return 0;
}

#ifdef _USE_CUDA
int reg_scale_amplitude_gpu(nifti_image *niftiImage, float min_value, float max_value)
{
//	reg_intensityRescale(niftiImage, min_value, max_value, -FLT_MAX, FLT_MAX);
	return 1;
}
#endif



int reg_gradient_jacobian_determinant(nifti_image *nodesGradientImage, nifti_image *controlPointImage, int image_size[], float image_spacing[], float weight)
{
	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);
	reg_bspline_jacobianDeterminantGradient<float>( controlPointImage,
								targetImage,
								nodesGradientImage,
								weight,
								1);
	free_nifti_image_except_data(targetImage);
	return 0;
}

#ifdef _USE_CUDA
int reg_gradient_jacobian_determinant_gpu(nifti_image *nodesGradientImage, nifti_image *controlPointImage, int image_size[], float image_spacing[], float weight)
{
//	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);
//	reg_bspline_jacobianDeterminantGradient<float>( controlPointImage,
//								targetImage,
//								nodesGradientImage,
//								weight,
//								1);
//	free_nifti_image_except_data(targetImage);
	return 1;
}
#endif


int reg_gradient_bending_energy(nifti_image *nodesGradientImage, nifti_image *controlPointImage, int image_size[], float image_spacing[], float weight)
{
	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);

//fprintf(stderr, "Target delta:       %f %f %f \n",targetImage->dx,targetImage->dy,targetImage->dz);
//fprintf(stderr, "Target size:        %d %d %d (nvox: %d) \n",targetImage->nx,targetImage->ny,targetImage->nz,targetImage->nvox);

//fprintf(stderr, "CP delta:           %f %f %f \n",controlPointImage->dx,controlPointImage->dy,controlPointImage->dz);
//fprintf(stderr, "CP size:            %d %d %d (nvox: %d) \n",controlPointImage->nx,controlPointImage->ny,controlPointImage->nz,controlPointImage->nvox);
//fprintf(stderr, "Weight:             %f \n",weight);

	reg_bspline_bendingEnergyGradient<float>(		controlPointImage,
								targetImage,
								nodesGradientImage,
								weight);
	free_nifti_image_except_data(targetImage);
	return 0;
}

#ifdef _USE_CUDA
int reg_gradient_bending_energy_gpu(nifti_image *nodesGradientImage, nifti_image *controlPointImage, int image_size[], float image_spacing[], float weight)
{
//	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);
//	reg_bspline_bendingEnergyGradient<float>(		controlPointImage,
//								targetImage,
//								nodesGradientImage,
//								weight);
//	free_nifti_image_except_data(targetImage);
	return 1;
}
#endif




nifti_image* reg_initialize_control_points(int control_point_size[], float gridSpacing[])
{
//fprintf(stderr,"Grid spacing: %f %f %f \n",gridSpacing[0],gridSpacing[1],gridSpacing[2]);
	int dim_cpp[8];
	dim_cpp[0]=5;
	dim_cpp[1]=control_point_size[0];
	dim_cpp[2]=control_point_size[1];
	dim_cpp[3]=control_point_size[2];
	dim_cpp[5]=3;
	dim_cpp[4]=dim_cpp[6]=dim_cpp[7]=1;
	nifti_image* controlPointImage = nifti_make_new_nim(dim_cpp, NIFTI_TYPE_FLOAT32, true);
	controlPointImage->cal_min=0;
	controlPointImage->cal_max=0;
	controlPointImage->pixdim[0]=1.0f;
	controlPointImage->pixdim[1]=controlPointImage->dx=gridSpacing[0];
	controlPointImage->pixdim[2]=controlPointImage->dy=gridSpacing[1];
	controlPointImage->pixdim[3]=controlPointImage->dz=gridSpacing[2];
	controlPointImage->pixdim[4]=controlPointImage->dt=1.0f;
	controlPointImage->pixdim[5]=controlPointImage->du=1.0f;
	controlPointImage->pixdim[6]=controlPointImage->dv=1.0f;
	controlPointImage->pixdim[7]=controlPointImage->dw=1.0f;
	controlPointImage->qform_code=1;
        controlPointImage->quatern_b=0.f;
        controlPointImage->quatern_c=0.f;
        controlPointImage->quatern_d=0.f;
        controlPointImage->qfac=1.f;
	controlPointImage->qoffset_x = -gridSpacing[0];
	controlPointImage->qoffset_y = -gridSpacing[1];
        controlPointImage->qoffset_z = -gridSpacing[2];
        controlPointImage->qto_xyz = nifti_quatern_to_mat44(0.f, 0.f, 0.f, -gridSpacing[0], -gridSpacing[1], -gridSpacing[2], gridSpacing[0], gridSpacing[1], gridSpacing[2], 1.f);
        controlPointImage->qto_ijk = nifti_mat44_inverse(controlPointImage->qto_xyz);
	controlPointImage->sform_code=0;
	return controlPointImage;
}


nifti_image *reg_initialize_image(int image_size[], float image_spacing[])
{
	int dim[8];
	dim[0]    = 3;
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *image = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
//	image->cal_min=0;
//	image->cal_max=0;
	image->pixdim[0]=1.0f;
	image->pixdim[1]=image->dx=image_spacing[0];
	image->pixdim[2]=image->dy=image_spacing[1];
	image->pixdim[3]=image->dz=image_spacing[2];
	image->pixdim[4]=image->dt=1.0f;
	image->pixdim[5]=image->du=1.0f;
	image->pixdim[6]=image->dv=1.0f;
	image->pixdim[7]=image->dw=1.0f;
	image->qform_code=1;
    image->quatern_b=0.f;
    image->quatern_c=0.f;
    image->quatern_d=0.f;
    image->qfac=1.f;
	image->qoffset_x = 0.f;
	image->qoffset_y = 0.f;
    image->qoffset_z = 0.f;

    image->qto_xyz = nifti_quatern_to_mat44(0.f, 0.f, 0.f,
	0, 0, 0,
	image_spacing[0], image_spacing[1], image_spacing[2], 1.f);
//fprintf(stderr,"\nTransformation matrix1: %f %f %f %f",image->qto_xyz.m[0][0],image->qto_xyz.m[0][1],image->qto_xyz.m[0][2],image->qto_xyz.m[0][3]);
//fprintf(stderr,"\nTransformation matrix2: %f %f %f %f",image->qto_xyz.m[1][0],image->qto_xyz.m[1][1],image->qto_xyz.m[1][2],image->qto_xyz.m[1][3]);
//fprintf(stderr,"\nTransformation matrix3: %f %f %f %f",image->qto_xyz.m[2][0],image->qto_xyz.m[2][1],image->qto_xyz.m[2][2],image->qto_xyz.m[2][3]);
//fprintf(stderr,"\nTransformation matrix4: %f %f %f %f",image->qto_xyz.m[3][0],image->qto_xyz.m[3][1],image->qto_xyz.m[3][2],image->qto_xyz.m[3][3]);
    image->qto_ijk = nifti_mat44_inverse(image->qto_xyz);
	image->sform_code=0;
	return image;	
}


nifti_image *reg_initialize_deformation_field(int image_size[], float image_spacing[])
{
	nifti_image *positionFieldImage = reg_initialize_image(image_size, image_spacing);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=3;
    positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    positionFieldImage->nbyper = sizeof(float);
	return positionFieldImage;
}


nifti_image *reg_initialize_image_4D(int image_size[], float image_spacing[])
{
	nifti_image *image = reg_initialize_image(image_size, image_spacing);
    image->dim[0]=image->ndim=5;
    image->dim[4]=image->nt=1;image->pixdim[4]=image->dt=1.0;
	image->dim[5]=image->nu=image_size[3];
    image->pixdim[5]=image->du=1.0;
    image->dim[6]=image->nv=1;image->pixdim[6]=image->dv=1.0;
    image->dim[7]=image->nw=1;image->pixdim[7]=image->dw=1.0;
    image->nvox=image->nx*image->ny*image->nz*image->nt*image->nu;
    image->datatype = NIFTI_TYPE_FLOAT32;
    image->nbyper = sizeof(float);
	return image;
}

int free_nifti_image_except_data(nifti_image *image)
{
	if( image->fname != NULL ) free(image->fname) ;
	if( image->iname != NULL ) free(image->iname) ;
	(void)nifti_free_extensions( image ) ;
	free(image) ;
	return 0;
}











#ifdef _USE_CUDA
unsigned int REG_resample_image_rigid_gpu(nifti_image *resampledImage, nifti_image *inputImage, mat44 *affineTransformation)
{
	cudaArray *inputImageArray_d = NULL; 
	float4    *positionFieldImageArray_d = NULL; 
	int       *resamplingMask_d = NULL; 
	float     *resampledImageArray_d = NULL; 

	int       *resamplingMask = NULL; 
	int       activeVoxelNumber = 0; 

    /* Allocate and initialize GPU memory */ 
    // This is a semiclever way to help not forget to free memory: 
    // (call alloc_record_add(..) after each malloc, as follows. )
    alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);

    // Input image:
    if(cudaCommon_allocateArrayToDevice<float>(&inputImageArray_d, inputImage->dim)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)inputImageArray_d,ALLOCTYPE_CUDA_ARRAY);
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&inputImageArray_d,inputImage)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    } 

    // Resampled (output) image: 
    if(cudaMalloc((void **)&resampledImageArray_d, resampledImage->nvox*sizeof(float)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    }  
    alloc_record_add(memory_record,(void*)resampledImageArray_d,ALLOCTYPE_CUDA); 

    // Mask 
	resamplingMask = (int *) malloc(resampledImage->nvox*sizeof(int)); 
	alloc_record_add(memory_record,(void*)resamplingMask,ALLOCTYPE_HOST); 
	for(unsigned int i=0; i<inputImage->nvox; i++)
		resamplingMask[i]=i; 
	activeVoxelNumber=inputImage->nvox;
    if(cudaMalloc((void **)&resamplingMask_d, resampledImage->nvox*sizeof(int)) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    }  
    alloc_record_add(memory_record,(void*)resamplingMask_d,ALLOCTYPE_CUDA); 
	if (cudaMemcpy(resamplingMask_d, (void*) resamplingMask, resampledImage->nvox*sizeof(int), cudaMemcpyHostToDevice)) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
	}

    // Position field 
    if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, resampledImage->dim)) { 
        alloc_record_destroy(memory_record); 
        return niftyrec_error_allocgpu;
    } 
    alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA); 

	/* Compute the position field */
    reg_affine_positionField_gpu(   affineTransformation, 
                                    resampledImage, 
                                    &positionFieldImageArray_d);

	/* Resample the image */
	reg_resampleSourceImage_gpu(		
	                    resampledImage,
						inputImage,
						&resampledImageArray_d,
						&inputImageArray_d,
						&positionFieldImageArray_d,
						&resamplingMask_d,
						resampledImage->nvox,
						0.0,
						0);  //background 
						
	/* Transfer data back to the host */					
	if(cudaCommon_transferFromDeviceToNifti(resampledImage, &resampledImageArray_d)) 
	    return niftyrec_error_transfergpu;


    int dim[8];
    dim[0]    = 4; 
    dim[1]    = inputImage->dim[1]; 
    dim[2]    = inputImage->dim[2];        
    dim[3]    = inputImage->dim[3]; 
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
    
    if (cudaMemcpy(((void *)tfieldImage->data), positionFieldImageArray_d, inputImage->nvox*sizeof(float)*4, cudaMemcpyDeviceToHost) != cudaSuccess) {
        alloc_record_destroy(memory_record); 
        return niftyrec_error_transfergpu;
    }

    nifti_set_filenames(tfieldImage,"/Users/spedemon/Desktop/tfield1.nii",0,0); 
    nifti_image_write(tfieldImage); 
    nifti_set_filenames(resampledImage,"/Users/spedemon/Desktop/in1.nii",0,0); 
    nifti_image_write(resampledImage); 



    /* Free */
	if (alloc_record_destroy(memory_record) != 0)
	    return STATUS_MEMORY_ERROR; 
	return STATUS_SUCCESS;  
}
#endif


unsigned int REG_resample_image_rigid_cpu(nifti_image *resampledImage, nifti_image *inputImage, mat44 *m)
{
    return 0; 
}






#ifdef _USE_CUDA
unsigned int REG_d_intensity_d_space_rigid_gpu(nifti_image *resultGradientImage, nifti_image *sourceImage, mat44 *affineTransformation) 
{
	int status;

	cudaArray *sourceImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	float4 *resultGradientArray_d=NULL;
	int *targetMask_d=NULL;

	int *targetMask;
	targetMask = (int *)malloc(sourceImage->nvox*sizeof(int));
	int activeVoxelNumber=0;
	for(unsigned int i=0; i<sourceImage->nvox; i++)
		targetMask[i]=i;
	activeVoxelNumber=sourceImage->nvox;

	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;

	// Index of the active voxel is stored
	int *targetMask_h; CUDA_SAFE_CALL(cudaMallocHost((void **)&targetMask_h, activeVoxelNumber*sizeof(int)));
	int *targetMask_h_ptr = &targetMask_h[0];
	for(unsigned int i=0;i<sourceImage->nvox;i++)
		{
		if(targetMask[i]!=-1) 
			*targetMask_h_ptr++=i;
		}
	CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, targetMask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&resultGradientArray_d, activeVoxelNumber*sizeof(float4)));

	// Compute the position field //
    reg_affine_positionField_gpu(       affineTransformation, 
                                        sourceImage, 
                                        &positionFieldImageArray_d); 

    // Compute gradient //
	reg_getSourceImageGradient_gpu(		
	                                    sourceImage,
						                sourceImage,
						                &sourceImageArray_d,
						                &positionFieldImageArray_d,
						                &resultGradientArray_d,
						                activeVoxelNumber,
						                0);

	if(cudaCommon_transferFromDeviceToNifti(resultGradientImage, &resultGradientArray_d)) return 1;

	free(targetMask);
	CUDA_SAFE_CALL(cudaFreeHost(targetMask_h));
	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&resultGradientArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	cudaCommon_free( (void **)&targetMask_d );

	status = STATUS_SUCCESS;
	return status;
}
#endif



unsigned int REG_d_intensity_d_space_rigid_cpu(nifti_image *gradientImage, nifti_image *inputImage, mat44 *m)
{
    return 0;
}




#ifdef _USE_CUDA
unsigned int REG_d_intensity_d_transformation_rigid_gpu(nifti_image *resultGradientImage, nifti_image *sourceImage, mat44 *affineTransformation) 
{
	int status;

	cudaArray *sourceImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	float4 *resultGradientArray_d=NULL;
	int *targetMask_d=NULL;

	int *targetMask;
	targetMask = (int *)malloc(sourceImage->nvox*sizeof(int));
	int activeVoxelNumber=0;
	for(unsigned int i=0; i<sourceImage->nvox; i++)
		targetMask[i]=i;
	activeVoxelNumber=sourceImage->nvox;

	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;

	// Index of the active voxel is stored
	int *targetMask_h; CUDA_SAFE_CALL(cudaMallocHost((void **)&targetMask_h, activeVoxelNumber*sizeof(int)));
	int *targetMask_h_ptr = &targetMask_h[0];
	for(unsigned int i=0;i<sourceImage->nvox;i++)
		{
		if(targetMask[i]!=-1) 
			*targetMask_h_ptr++=i;
		}
	CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, targetMask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&resultGradientArray_d, activeVoxelNumber*sizeof(float4)));

	// Compute the position field //
    reg_affine_positionField_gpu(       affineTransformation, 
                                        sourceImage, 
                                        &positionFieldImageArray_d); 

    // Compute gradient //
	reg_getSourceImageGradient_gpu(		
	                                    sourceImage,
						                sourceImage,
						                &sourceImageArray_d,
						                &positionFieldImageArray_d,
						                &resultGradientArray_d,
						                activeVoxelNumber,
						                0);

	if(cudaCommon_transferFromDeviceToNifti(resultGradientImage, &resultGradientArray_d)) return 1;

	free(targetMask);
	CUDA_SAFE_CALL(cudaFreeHost(targetMask_h));
	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&resultGradientArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	cudaCommon_free( (void **)&targetMask_d );

	status = STATUS_SUCCESS;
	return status;
}
#endif



unsigned int REG_d_intensity_d_transformation_rigid_cpu(nifti_image *gradientImage, nifti_image *inputImage, mat44 *m)
{
    return 0;
}


