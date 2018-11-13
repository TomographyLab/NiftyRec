/*
 *  _reg_cudaCommon.cpp
 *  
 *
 *  Created by Marc Modat on 25/03/2009.
 *  Copyright (c) 2009, University College London. All rights reserved.
 *  Centre for Medical Image Computing (CMIC)
 *  See the LICENSE.txt file in the nifty_reg root folder
 *
 */
#ifdef _USE_CUDA

#ifndef _REG_CUDACOMMON_CPP
#define _REG_CUDACOMMON_CPP

#include "_reg_cudaCommon.h"




/* ******************************** */
/* ******************************** */
template <class DTYPE, class NIFTI_TYPE>
int cudaCommon_transferNiftiToArrayOnDevice1(DTYPE **array_d, nifti_image *img)
{
	if(sizeof(DTYPE)!=sizeof(NIFTI_TYPE)){
		fprintf(stderr, "ERROR:\tcudaCommon_transferNiftiToArrayOnDevice:\n");
		fprintf(stderr, "ERROR:\tThe host and device arrays are of different types.\n");
		return 1;
	}
	else{
		const unsigned int memSize = img->dim[1] * img->dim[2] * img->dim[3] * sizeof(DTYPE);
		NIFTI_TYPE *array_h=static_cast<NIFTI_TYPE *>(img->data);
		CUDA_SAFE_CALL(cudaMemcpy(*array_d, array_h, memSize, cudaMemcpyHostToDevice));
	}
	return 0;
}
/* ******************************** */
template <class DTYPE>
int cudaCommon_transferNiftiToArrayOnDevice(DTYPE **array_d, nifti_image *img)
{
	if( sizeof(DTYPE)==sizeof(float4) ){
		if( (img->datatype!=NIFTI_TYPE_FLOAT32) || (img->dim[5]<2) || (img->dim[4]>1)){
			fprintf(stderr, "ERROR:\tcudaCommon_transferNiftiToDevice:\n");
			fprintf(stderr, "ERROR:\tThe specified image is not a single precision deformation field image\n");
			return 1;
		}
		float *niftiImgValues = static_cast<float *>(img->data);
		float4 *array_h=(float4 *)calloc(img->nx*img->ny*img->nz,sizeof(float4));
		const int voxelNumber = img->nx*img->ny*img->nz;
		for(int i=0; i<voxelNumber; i++)
			array_h[i].x= *niftiImgValues++;
		if(img->dim[5]>=2){
			for(int i=0; i<voxelNumber; i++)
				array_h[i].y= *niftiImgValues++;
		}
		if(img->dim[5]>=3){
			for(int i=0; i<voxelNumber; i++)
				array_h[i].z= *niftiImgValues++;
		}
		if(img->dim[5]>=4){
			for(int i=0; i<voxelNumber; i++)
				array_h[i].w= *niftiImgValues++;
		}
		CUDA_SAFE_CALL(cudaMemcpy(*array_d, array_h, img->nx*img->ny*img->nz*sizeof(float4), cudaMemcpyHostToDevice));
		free(array_h);
	}
	else{ // All these else could be removed but the nvcc compiler would warn for unreachable statement
		switch(img->datatype){
			case NIFTI_TYPE_UINT8:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,unsigned char>(array_d, img);
			case NIFTI_TYPE_INT8:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,char>(array_d, img);
			case NIFTI_TYPE_UINT16:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,unsigned short>(array_d, img);
			case NIFTI_TYPE_INT16:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,short>(array_d, img);
			case NIFTI_TYPE_UINT32:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,unsigned int>(array_d, img);
			case NIFTI_TYPE_INT32:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,int>(array_d, img);
			case NIFTI_TYPE_FLOAT32:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,float>(array_d, img);
			case NIFTI_TYPE_FLOAT64:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,double>(array_d, img);
			default:
				fprintf(stderr, "ERROR:\tcudaCommon_transferNiftiToArrayOnDevice:\n");
				fprintf(stderr, "ERROR:\tThe image data type is not supported\n");
				return 1;
		}
	}
	return 0;
}
template int cudaCommon_transferNiftiToArrayOnDevice<char>(char **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<unsigned char>(unsigned char **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<short>(short **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<unsigned short>(unsigned short **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<int>(int **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<unsigned int>(unsigned int **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<float>(float **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<double>(double **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<float4>(float4 **, nifti_image *); // for deformation field
/* ******************************** */
/* ******************************** */
template <class DTYPE, class NIFTI_TYPE>
int cudaCommon_transferNiftiToArrayOnDevice1(cudaArray **cuArray_d, nifti_image *img)
{
	if(sizeof(DTYPE)!=sizeof(NIFTI_TYPE)){
		fprintf(stderr, "ERROR:\tcudaCommon_transferNiftiToArrayOnDevice:\n");
		fprintf(stderr, "ERROR:\tThe host and device arrays are of different types.\n");
		return 1;
	}
	else{
		const cudaExtent volumeSize = make_cudaExtent(img->dim[1], img->dim[2], img->dim[3]);
		NIFTI_TYPE *array_h=static_cast<NIFTI_TYPE *>(img->data);

		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr = make_cudaPitchedPtr((void *) array_h, volumeSize.width*sizeof(DTYPE), volumeSize.width, volumeSize.height);
		copyParams.dstArray = *cuArray_d;
		copyParams.extent = volumeSize;
		copyParams.kind = cudaMemcpyHostToDevice;
		CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
		
	}
	return 0;
}
/* ******************************** */
template <class DTYPE>
int cudaCommon_transferNiftiToArrayOnDevice(cudaArray **cuArray_d, nifti_image *img)
{
	if( sizeof(DTYPE)==sizeof(float4) ){
		if( (img->datatype!=NIFTI_TYPE_FLOAT32) || (img->dim[5]<2) || (img->dim[4]>1) ){
			fprintf(stderr, "ERROR:\tcudaCommon_transferNiftiToDevice:\n");
			fprintf(stderr, "ERROR:\tThe specified image is not a single precision deformation field image\n");
			return 1;
		}
		float *niftiImgValues = static_cast<float *>(img->data);
		float4 *array_h=(float4 *)calloc(img->nx*img->ny*img->nz,sizeof(float4));
		for(int i=0; i<img->nx*img->ny*img->nz; i++)
			(array_h++)->x= *niftiImgValues++;
		if(img->dim[5]>=2){
			for(int i=0; i<img->nx*img->ny*img->nz; i++)
				(array_h++)->y= *niftiImgValues++;
		}
		if(img->dim[5]>=3){
			for(int i=0; i<img->nx*img->ny*img->nz; i++)
				(array_h++)->z= *niftiImgValues++;
		}
		if(img->dim[5]==3){
			for(int i=0; i<img->nx*img->ny*img->nz; i++)
				(array_h++)->w= *niftiImgValues++;
		}
		const cudaExtent volumeSize = make_cudaExtent(img->dim[1], img->dim[2], img->dim[3]);
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr = make_cudaPitchedPtr((void *) array_h, volumeSize.width*sizeof(DTYPE), volumeSize.width, volumeSize.height);
		copyParams.dstArray = *cuArray_d;
		copyParams.extent = volumeSize;
		copyParams.kind = cudaMemcpyHostToDevice;
		CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
		free(array_h);
	}
	else{ // All these else could be removed but the nvcc compiler would warn for unreachable statement 
		switch(img->datatype){
			case NIFTI_TYPE_UINT8:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,unsigned char>(cuArray_d, img);
			case NIFTI_TYPE_INT8:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,char>(cuArray_d, img);
			case NIFTI_TYPE_UINT16:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,unsigned short>(cuArray_d, img);
			case NIFTI_TYPE_INT16:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,short>(cuArray_d, img);
			case NIFTI_TYPE_UINT32:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,unsigned int>(cuArray_d, img);
			case NIFTI_TYPE_INT32:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,int>(cuArray_d, img);
			case NIFTI_TYPE_FLOAT32:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,float>(cuArray_d, img);
			case NIFTI_TYPE_FLOAT64:
				return cudaCommon_transferNiftiToArrayOnDevice1<DTYPE,double>(cuArray_d, img);
			default:
				fprintf(stderr, "ERROR:\tcudaCommon_transferNiftiToArrayOnDevice:\n");
				fprintf(stderr, "ERROR:\tThe image data type is not supported\n");
				return 1;
		}
	}
	return 0;
}
template int cudaCommon_transferNiftiToArrayOnDevice<char>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<unsigned char>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<short>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<unsigned short>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<int>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<unsigned int>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<float>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<double>(cudaArray **, nifti_image *);
template int cudaCommon_transferNiftiToArrayOnDevice<float4>(cudaArray **, nifti_image *); // for deformation field
/* ******************************** */
/* ******************************** */
template <class DTYPE>
int cudaCommon_allocateArrayToDevice(cudaArray **cuArray_d, int *dim)
{
	const cudaExtent volumeSize = make_cudaExtent(dim[1], dim[2], dim[3]);
	cudaChannelFormatDesc texDesc = cudaCreateChannelDesc<DTYPE>();
	CUDA_SAFE_CALL(cudaMalloc3DArray(cuArray_d, &texDesc, volumeSize));
	return 0;
}
template int cudaCommon_allocateArrayToDevice<char>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<unsigned char>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<short>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<unsigned short>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<int>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<unsigned int>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<float>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<double>(cudaArray **, int *);
template int cudaCommon_allocateArrayToDevice<float4>(cudaArray **, int *); // for deformation field
/* ******************************** */
/* ******************************** */
template <class DTYPE>
int cudaCommon_allocateArrayToDevice(DTYPE **array_d, int *dim)
{
	const unsigned int memSize = dim[1] * dim[2] * dim[3] * sizeof(DTYPE);
	CUDA_SAFE_CALL(cudaMalloc((void **)array_d, memSize));
	return 0;
}
template int cudaCommon_allocateArrayToDevice<char>(char **, int *);
template int cudaCommon_allocateArrayToDevice<unsigned char>(unsigned char **, int *);
template int cudaCommon_allocateArrayToDevice<short>(short **, int *);
template int cudaCommon_allocateArrayToDevice<unsigned short>(unsigned short **, int *);
template int cudaCommon_allocateArrayToDevice<int>(int **, int *);
template int cudaCommon_allocateArrayToDevice<unsigned int>(unsigned int **, int *);
template int cudaCommon_allocateArrayToDevice<float>(float **, int *);
template int cudaCommon_allocateArrayToDevice<double>(double **, int *);
template int cudaCommon_allocateArrayToDevice<float4>(float4 **, int *); // for deformation field
/* ******************************** */
/* ******************************** */
template <class DTYPE, class NIFTI_TYPE>
int cudaCommon_transferFromDeviceToNifti1(nifti_image *img, DTYPE **array_d)
{
	if(sizeof(DTYPE)!=sizeof(NIFTI_TYPE)){
		fprintf(stderr, "ERROR:\tcudaCommon_transferFromDeviceToNifti:\n");
		fprintf(stderr, "ERROR:\tThe host and device arrays are of different types.\n");
		return 1;
	}
	else{
		NIFTI_TYPE *array_h=static_cast<NIFTI_TYPE *>(img->data);
		CUDA_SAFE_CALL(cudaMemcpy(array_h, *array_d, img->nvox*sizeof(DTYPE), cudaMemcpyDeviceToHost));
	}
	return 0;
}
/* ******************************** */
template <class DTYPE>
int cudaCommon_transferFromDeviceToNifti(nifti_image *img, DTYPE **array_d)
{
	if(sizeof(DTYPE)==sizeof(float4)){
        // A nifti 5D volume is expected
		if(img->dim[0]<5 || img->dim[4]>1 || img->dim[5]<2 || img->datatype!=NIFTI_TYPE_FLOAT32){
			fprintf(stderr, "ERROR:\tcudaCommon_transferFromDeviceToNifti:\n");
			fprintf(stderr, "ERROR:\tThe nifti image is not a 5D volume.\n");
			return 1;
		}
		const int voxelNumber = img->nx*img->ny*img->nz;
		float4 *array_h = (float4 *)calloc(voxelNumber, sizeof(float4));
		CUDA_SAFE_CALL(cudaMemcpy(array_h, *array_d, voxelNumber*sizeof(float4), cudaMemcpyDeviceToHost));
		float *niftiImgValues = static_cast<float *>(img->data);
		for(int i=0; i<voxelNumber; i++)
			*niftiImgValues++ = array_h[i].x;
		if(img->dim[5]>=2){
			for(int i=0; i<voxelNumber; i++)
				*niftiImgValues++ = array_h[i].y;
		}
		if(img->dim[5]>=3){
			for(int i=0; i<voxelNumber; i++)
				*niftiImgValues++ = array_h[i].z;
		}
		if(img->dim[5]>=4){
			for(int i=0; i<voxelNumber; i++)
				*niftiImgValues++ = array_h[i].w;
		}
		free(array_h);
		
		return 0;
	}
	else{
		switch(img->datatype){
			case NIFTI_TYPE_UINT8:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,unsigned char>(img, array_d);
			case NIFTI_TYPE_INT8:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,char>(img, array_d);
			case NIFTI_TYPE_UINT16:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,unsigned short>(img, array_d);
			case NIFTI_TYPE_INT16:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,short>(img, array_d);
			case NIFTI_TYPE_UINT32:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,unsigned int>(img, array_d);
			case NIFTI_TYPE_INT32:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,int>(img, array_d);
			case NIFTI_TYPE_FLOAT32:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,float>(img, array_d);
			case NIFTI_TYPE_FLOAT64:
				return cudaCommon_transferFromDeviceToNifti1<DTYPE,double>(img, array_d);
			default:
				fprintf(stderr, "ERROR:\tcudaCommon_transferFromDeviceToNifti:\n");
				fprintf(stderr, "ERROR:\tThe image data type is not supported\n");
				return 1;
		}
	}
}
template int cudaCommon_transferFromDeviceToNifti<char>(nifti_image *, char **);
template int cudaCommon_transferFromDeviceToNifti<unsigned char>(nifti_image *, unsigned char **);
template int cudaCommon_transferFromDeviceToNifti<short>(nifti_image *, short **);
template int cudaCommon_transferFromDeviceToNifti<unsigned short>(nifti_image *, unsigned short **);
template int cudaCommon_transferFromDeviceToNifti<int>(nifti_image *, int **);
template int cudaCommon_transferFromDeviceToNifti<unsigned int>(nifti_image *, unsigned int **);
template int cudaCommon_transferFromDeviceToNifti<float>(nifti_image *, float **);
template int cudaCommon_transferFromDeviceToNifti<double>(nifti_image *, double **);
template int cudaCommon_transferFromDeviceToNifti<float4>(nifti_image *, float4 **); // for deformation field
/* ******************************** */
/* ******************************** */
void cudaCommon_free(cudaArray **cuArray_d){
	cudaFreeArray(*cuArray_d);
	return;
}
/* ******************************** */
/* ******************************** */
void cudaCommon_free(void **array_d){
	cudaFree(*array_d);
	return;
}
/* ******************************** */
/* ******************************** */
#endif
#endif
