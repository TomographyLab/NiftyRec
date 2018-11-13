/*
 *  _reg_resampling_gpu.cu
 *  
 *
 *  Created by Marc Modat on 24/03/2009.
 *  Copyright (c) 2009, University College London. All rights reserved.
 *  Centre for Medical Image Computing (CMIC)
 *  See the LICENSE.txt file in the nifty_reg root folder
 *
 */

#ifndef _REG_RESAMPLING_GPU_CU
#define _REG_RESAMPLING_GPU_CU
#define REG_RESAMPLING_LINEAR 0
#define REG_RESAMPLING_POINT  1

#include "_reg_resampling_gpu.h"
#include "_reg_resampling_gpu_kernels.cu"

void reg_resampleSourceImage_gpu(       nifti_image *resultImage,
                                        nifti_image *sourceImage,
                                        float **resultImageArray_d,
                                        cudaArray **sourceImageArray_d,
                                        float4 **positionFieldImageArray_d,
                                        int **mask_d,
                                        int activeVoxelNumber,
                                        float sourceBGValue, 
                                        unsigned int resampling_mode)
{
        int3 sourceDim = make_int3(sourceImage->nx, sourceImage->ny, sourceImage->nz);

        CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_SourceDim,&sourceDim,sizeof(int3)));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_PaddingValue,&sourceBGValue,sizeof(float)));
        CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ActiveVoxelNumber,&activeVoxelNumber,sizeof(int)));

        //Bind source image array to a 3D texture
        sourceTexture.normalized = true; 

        if (resampling_mode == REG_RESAMPLING_LINEAR)
            sourceTexture.filterMode = cudaFilterModeLinear; 
        if (resampling_mode == REG_RESAMPLING_POINT)   
            sourceTexture.filterMode = cudaFilterModePoint; 
            
        sourceTexture.addressMode[0] = cudaAddressModeWrap;
        sourceTexture.addressMode[1] = cudaAddressModeWrap;
        sourceTexture.addressMode[2] = cudaAddressModeWrap;

        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        CUDA_SAFE_CALL(cudaBindTextureToArray(sourceTexture, *sourceImageArray_d, channelDesc));

        //Bind positionField to texture
        CUDA_SAFE_CALL(cudaBindTexture(0, positionFieldTexture, *positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));

        //Bind positionField to texture
        CUDA_SAFE_CALL(cudaBindTexture(0, maskTexture, *mask_d, activeVoxelNumber*sizeof(int)));

        // Bind the real to voxel matrix to texture
        mat44 *sourceMatrix;
        if(sourceImage->sform_code>0)
                sourceMatrix=&(sourceImage->sto_ijk);
        else sourceMatrix=&(sourceImage->qto_ijk);
        float4 *sourceRealToVoxel_h;CUDA_SAFE_CALL(cudaMallocHost((void **)&sourceRealToVoxel_h, 3*sizeof(float4)));
        float4 *sourceRealToVoxel_d;
        CUDA_SAFE_CALL(cudaMalloc((void **)&sourceRealToVoxel_d, 3*sizeof(float4)));
        for(int i=0; i<3; i++){
                sourceRealToVoxel_h[i].x=sourceMatrix->m[i][0];
                sourceRealToVoxel_h[i].y=sourceMatrix->m[i][1];
                sourceRealToVoxel_h[i].z=sourceMatrix->m[i][2];
                sourceRealToVoxel_h[i].w=sourceMatrix->m[i][3];
        }
        

//        fprintf(stderr,"@@@@@@@@  Resampling matrix source->sto_ijk: @@@@@@@@@\n");
//        for (int i=0; i<4; i++)
//			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",sourceMatrix->m[i][0],sourceMatrix->m[i][1],sourceMatrix->m[i][2],sourceMatrix->m[i][3]);
//        fprintf(stderr,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"); 
            
            
        
        CUDA_SAFE_CALL(cudaMemcpy(sourceRealToVoxel_d, sourceRealToVoxel_h, 3*sizeof(float4), cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaFreeHost((void *)sourceRealToVoxel_h));
        CUDA_SAFE_CALL(cudaBindTexture(0, sourceMatrixTexture, sourceRealToVoxel_d, 3*sizeof(float4)));

        const unsigned int Grid_reg_resampleSourceImage = (unsigned int)ceil((float)activeVoxelNumber/(float)Block_reg_resampleSourceImage);
        dim3 B1(Block_reg_resampleSourceImage,1,1);
        dim3 G1(Grid_reg_resampleSourceImage,1,1);

        reg_resampleSourceImage_kernel <<< G1, B1 >>> (*resultImageArray_d);
        CUDA_SAFE_CALL(cudaThreadSynchronize());
#if _VERBOSE
        printf("[VERBOSE] reg_resampleSourceImage_kernel kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
                cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
        // Unbind 
        CUDA_SAFE_CALL(cudaUnbindTexture(sourceTexture));
        CUDA_SAFE_CALL(cudaUnbindTexture(positionFieldTexture));
        CUDA_SAFE_CALL(cudaUnbindTexture(maskTexture));
        CUDA_SAFE_CALL(cudaUnbindTexture(sourceMatrixTexture));

        cudaFree(sourceRealToVoxel_d);
}








void reg_getSourceImageGradient_gpu(	nifti_image *targetImage,
					nifti_image *sourceImage,
					cudaArray **sourceImageArray_d,
					float4 **positionFieldImageArray_d,
					float4 **resultGradientArray_d,
                    int activeVoxelNumber, 
                    unsigned int resampling_mode )
{
	int3 sourceDim = make_int3(sourceImage->nx, sourceImage->ny, sourceImage->nz);

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_SourceDim, &sourceDim, sizeof(int3)));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ActiveVoxelNumber, &activeVoxelNumber, sizeof(int)));

	//Bind source image array to a 3D texture
	sourceTexture.normalized = true;
    if (resampling_mode == REG_RESAMPLING_LINEAR)
        sourceTexture.filterMode = cudaFilterModeLinear; 
    if (resampling_mode == REG_RESAMPLING_POINT)   
        sourceTexture.filterMode = cudaFilterModePoint; 
	sourceTexture.addressMode[0] = cudaAddressModeWrap;
	sourceTexture.addressMode[1] = cudaAddressModeWrap;
	sourceTexture.addressMode[2] = cudaAddressModeWrap;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	CUDA_SAFE_CALL(cudaBindTextureToArray(sourceTexture, *sourceImageArray_d, channelDesc));

	//Bind positionField to texture
	cudaBindTexture(0, positionFieldTexture, *positionFieldImageArray_d, activeVoxelNumber*sizeof(float4));

	// Bind the real to voxel matrix to texture
	mat44 *sourceMatrix;
	if(sourceImage->sform_code>0)
		sourceMatrix=&(sourceImage->sto_ijk);
	else sourceMatrix=&(sourceImage->qto_ijk);
	float4 *sourceRealToVoxel_h;CUDA_SAFE_CALL(cudaMallocHost((void **)&sourceRealToVoxel_h, 3*sizeof(float4)));
	float4 *sourceRealToVoxel_d;
	CUDA_SAFE_CALL(cudaMalloc((void **)&sourceRealToVoxel_d, 3*sizeof(float4)));
	for(int i=0; i<3; i++){
		sourceRealToVoxel_h[i].x=sourceMatrix->m[i][0];
		sourceRealToVoxel_h[i].y=sourceMatrix->m[i][1];
		sourceRealToVoxel_h[i].z=sourceMatrix->m[i][2];
		sourceRealToVoxel_h[i].w=sourceMatrix->m[i][3];
	}
	CUDA_SAFE_CALL(cudaMemcpy(sourceRealToVoxel_d, sourceRealToVoxel_h, 3*sizeof(float4), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaFreeHost((void *)sourceRealToVoxel_h));
	CUDA_SAFE_CALL(cudaBindTexture(0, sourceMatrixTexture, sourceRealToVoxel_d, 3*sizeof(float4)));

	const unsigned int Grid_reg_getSourceImageGradient = (unsigned int)ceil((float)activeVoxelNumber/(float)Block_reg_getSourceImageGradient);
	dim3 B1(Block_reg_getSourceImageGradient,1,1);
	dim3 G1(Grid_reg_getSourceImageGradient,1,1);

	reg_getSourceImageGradient_kernel <<< G1, B1 >>> (*resultGradientArray_d);
	CUDA_SAFE_CALL(cudaThreadSynchronize());
#if _VERBOSE
	printf("[VERBOSE] reg_getSourceImageGradient kernel: %s - Grid size [%i %i %i] - Block size [%i %i %i]\n",
	       cudaGetErrorString(cudaGetLastError()),G1.x,G1.y,G1.z,B1.x,B1.y,B1.z);
#endif
	cudaFree(sourceRealToVoxel_d);
}

#endif
