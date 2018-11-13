/*
 *  _et_convolve2D.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_convolve2D.h"

int convolve2D_image(float *image_data, int image_size_x, int image_size_y, float *kernel_data, int kernel_size_x, int kernel_size_y, float *return_data, float background)
{
    return 0;
}

int et_convolve2D(nifti_image *inputImage, nifti_image *kernelImage, nifti_image *outImage, float background)
{
    int status = 1; 
    const int dataH = inputImage->nx; 
    const int dataW = inputImage->ny; 
    const int kernelH = kernelImage->nx; 
    const int kernelW = kernelImage->ny; 

    const int kernelX = (kernelH-1)/2;
    const int kernelY = (kernelW-1)/2;

    const int n_slices = inputImage->nz; 
    const int data_slice_size = dataH * dataW;
    const int kernel_slice_size = kernelH * kernelW;

    float *ptr_Data, *ptr_Kernel, ptr_Result; 

    float *convolved_slice = (float*) malloc(dataH*dataW*sizeof(float));

    for (int slice=0; slice<n_slices; slice++) 
        {
        ptr_Data = (((float*)inputImage->data) + slice * data_slice_size);
        ptr_Kernel = ((float*)kernelImage->data + slice * kernel_slice_size);
        convolve2D_image(ptr_Data, dataH, dataW, ptr_Kernel, kernelH, kernelW, convolved_slice, background);
        memcpy(ptr_Data, convolved_slice, dataH*dataW*sizeof(float));
        }
    free(convolved_slice);
    return status;
}
