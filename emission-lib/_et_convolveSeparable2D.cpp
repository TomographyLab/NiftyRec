/*
 *  _et_convolveSeparable2D.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_convolveSeparable2D.h"

extern "C" void convolutionRow(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR)
{
    for(int y = 0; y < imageH; y++)
        for(int x = 0; x < imageW; x++){
            float sum = 0;
            for(int k = -kernelR; k <= kernelR; k++){
                int d = x + k;
                if(d >= 0 && d < imageW)
                    sum += h_Src[y * imageW + d] * h_Kernel[kernelR - k];
            }
            h_Dst[y * imageW + x] = sum;
        }
}


extern "C" void convolutionColumn(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR)
{
    for(int y = 0; y < imageH; y++)
        for(int x = 0; x < imageW; x++){
            float sum = 0;
            for(int k = -kernelR; k <= kernelR; k++){
                int d = y + k;
                if(d >= 0 && d < imageH)
                    sum += h_Src[d * imageW + x] * h_Kernel[kernelR - k];
            }
            h_Dst[y * imageW + x] = sum;
        }
}


int et_convolveSeparable2D(nifti_image* inputImage, float *kernelSeparated, int kernelLengthH, int kernelLengthW, nifti_image *outputImage, float background)
{
    int imageH = inputImage->ny;
    int imageW = inputImage->nx;
    int n_slices = inputImage->nz;
    int kernelRadiusH = (kernelLengthH-1)/2;
    int kernelRadiusW = (kernelLengthW-1)/2;

    int image_slice_size = imageH * imageW;
    int kernel_slice_size = kernelLengthW + kernelLengthH;

    float *buffer = (float*) malloc(imageH*imageW*sizeof(float));

    float *input, *output, *kernel_row, *kernel_column;

    //Convolve slices one by one
    for (int slice=0; slice<n_slices; slice++)
        {
        // determine slice pointers
        input = (float*)inputImage->data + slice * image_slice_size;
        output = (float*)outputImage->data + slice * image_slice_size; 
        kernel_row = kernelSeparated + slice * kernel_slice_size;
        kernel_column = kernelSeparated + slice * kernel_slice_size + kernelLengthW;

        // convolve slice
        convolutionRow(buffer,input,kernel_row,imageW,imageH,kernelRadiusW);
        convolutionColumn(output,buffer,kernel_column,imageW,imageH,kernelRadiusH);
        }

    free(buffer);

    return 0;
}




/* FIXME Multi-threaded version is broken: 
#include "_et_convolveSeparable2D.h"
#include <pthread.h>
#define N_THREADS 4

typedef struct {
    float *input_data;
    float *buffer; 
    float *kernelSeparated;
    float *output_data; 
    int start_slice;
    int end_slice;
    int imageW;
    int imageH;
    int kernelRadiusW;
    int kernelRadiusH;
    int kernelLengthW;
    int kernelLengthH;
} thread_param;


extern "C" void convolutionRow(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR)
{
    for(int y = 0; y < imageH; y++)
        for(int x = 0; x < imageW; x++){
            float sum = 0;
            for(int k = -kernelR; k <= kernelR; k++){
                int d = x + k;
                if(d >= 0 && d < imageW)
                    sum += h_Src[y * imageW + d] * h_Kernel[kernelR - k];
            }
            h_Dst[y * imageW + x] = sum;
        }
}


extern "C" void convolutionColumn(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR)
{
    for(int y = 0; y < imageH; y++)
        for(int x = 0; x < imageW; x++){
            float sum = 0;
            for(int k = -kernelR; k <= kernelR; k++){
                int d = y + k;
                if(d >= 0 && d < imageH)
                    sum += h_Src[d * imageW + x] * h_Kernel[kernelR - k];
            }
            h_Dst[y * imageW + x] = sum;
        }
}


void *convolve_stack_of_slices(void *arg)
{
    thread_param *tp=(thread_param *)arg;
    float *input_data = tp->input_data;
    float *buffer = tp-> buffer;
    float *kernelSeparated = tp->kernelSeparated;
    float *output_data = tp->output_data;
    int start_slice = tp->start_slice;
    int end_slice = tp->end_slice;
    int imageW = tp->imageW;
    int imageH = tp->imageH;
    int kernelRadiusW = tp->kernelRadiusW;
    int kernelRadiusH = tp->kernelRadiusH;
    int kernelLengthW = tp->kernelLengthW;
    int kernelLengthH = tp->kernelLengthH;

    float *input, *output, *kernel_row, *kernel_column;
    int image_slice_size = imageH * imageW;
    int kernel_slice_size = kernelLengthW + kernelLengthH;

    for (int slice=start_slice; slice<end_slice; slice++) 
        {
        // determine slice pointers 
        input = (float*)input_data + slice * image_slice_size;
        output = (float*)output_data + slice * image_slice_size; 
        kernel_row = kernelSeparated + slice * kernel_slice_size;
        kernel_column = kernelSeparated + slice * kernel_slice_size + kernelLengthW;

        // convolve slice 
        convolutionRow(buffer,input,kernel_row,imageW,imageH,kernelRadiusW);
        convolutionColumn(output,buffer,kernel_column,imageW,imageH,kernelRadiusH);
        }
}


int et_convolveSeparable2D(nifti_image* inputImage, float *kernelSeparated, int kernelLengthH, int kernelLengthW, nifti_image *outputImage, float background)
{
    int threaded = 1;

    int imageH = inputImage->ny;
    int imageW = inputImage->nx;
    int n_slices = inputImage->nz;
    int kernelRadiusH = (kernelLengthH-1)/2;
    int kernelRadiusW = (kernelLengthW-1)/2;
    float *buffer = (float*) malloc(imageH*imageW*sizeof(float));

    //Convolve slices one by one
    pthread_t *threads=(pthread_t *)malloc(N_THREADS*sizeof(*threads)); 
    thread_param *tp=(thread_param *)malloc(sizeof(thread_param)*N_THREADS); 
    int n_slices_per_thread = ceil(n_slices/N_THREADS);

    for (int i=0; i<N_THREADS; i++) 
        {
        tp[i].input_data = (float*)inputImage->data;
        tp[i].buffer = buffer;
        tp[i].kernelSeparated = kernelSeparated;
        tp[i].output_data = (float*)outputImage->data;
        tp[i].start_slice = i*n_slices_per_thread; 
        tp[i].end_slice   = tp[i].start_slice+n_slices_per_thread; 
        if (tp[i].end_slice>=n_slices) tp[i].end_slice=n_slices; 
        tp[i].imageW = imageW; 
        tp[i].imageH = imageH; 
        tp[i].kernelRadiusW = kernelRadiusW; 
        tp[i].kernelRadiusH = kernelRadiusH; 
        tp[i].kernelLengthW = kernelLengthW; 
        tp[i].kernelLengthH = kernelLengthH; 
        if (threaded)
            pthread_create(&threads[i],NULL,convolve_stack_of_slices,(void *)(tp+i));
        else
            convolve_stack_of_slices((void *)(tp+i));
        }
        //wait for all threads 
        if (threaded)
            for (int i=0; i<N_THREADS; i++) 
                pthread_join(threads[i],NULL);
    free(threads);
    free(tp);  
    free(buffer);

    return 0;
}
*/

