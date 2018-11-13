/*
 *  _et_line_integral.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_line_integral.h"
//#include <pthread.h>
//#include <time.h>
#define N_THREADS 1

//long int timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p)
//{
//  return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
//           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
//}


typedef struct {
 int nx,ny,nz,start_x,end_x; 
 float *sino_data; 
 float *input_data; 
} thread_param;

void line_integral_thread(void *arg)
{
    thread_param *p=(thread_param *)arg;
    float *sino_data=p->sino_data;
    float *input_data=p->input_data;
    for(int x=p->start_x; x<p->end_x; x++)  
        for(int y=0; y<p->ny; y++)  
            for(int z=0; z<p->nz; z++) 
                sino_data[x + y*p->nx] = sino_data[x + y*p->nx] + input_data[x + y*p->nx + z*p->nx*p->ny];  
}

void et_line_integral(nifti_image *inputImage, nifti_image *sinoImage, int cam, float background)
{
//    struct timespec start, end;
//    clock_gettime(CLOCK_MONOTONIC, &start);
    int threaded=0;

    float *sino_data = (float *) (sinoImage->data) + cam*inputImage->nx*inputImage->ny ;
    float *input_data = (float *) (inputImage->data);
    memset((void*) sino_data,0,inputImage->nx*inputImage->ny*sizeof(float));
    int nx = inputImage->nx;
    int ny = inputImage->ny;
    int nz = inputImage->nz;
//    pthread_t *threads=(pthread_t *)malloc(N_THREADS*sizeof(*threads)); 
    thread_param *p=(thread_param *)malloc(sizeof(thread_param)*N_THREADS); 
    int n_lines_per_thread = ceil((double)nx/N_THREADS);
    for (int i=0; i<N_THREADS; i++) 
        {
        p[i].nx = nx; 
        p[i].ny = ny; 
        p[i].nz = nz; 
        p[i].start_x = i*n_lines_per_thread; 
        p[i].end_x   = p[i].start_x+n_lines_per_thread; 
        if (p[i].end_x>=nx) p[i].end_x=nx; 
        p[i].sino_data = sino_data;
        p[i].input_data = input_data;  
 //       if (threaded)
 //           pthread_create(&threads[i],NULL,line_integral_thread,(void *)(p+i));
 //       else
            line_integral_thread((void *)(p+i));
        }
        //wait for all threads 
 //       if (threaded)
 //           for (int i=0; i<N_THREADS; i++) 
 //               pthread_join(threads[i],NULL);
 //   free(threads);
    free(p);  

//    clock_gettime(CLOCK_MONOTONIC, &end);
//    long int timeElapsed = timespecDiff(&end, &start);
//    fprintf(stderr,"Time to compute the line integrals: %d usec. \n",timeElapsed/1000);
}    






