/*
 *  _et_line_backproject_attenuated.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et_line_backproject_attenuated.h"
#ifdef _OPENMP
#include "omp.h"
#endif

void et_line_backproject_attenuated(nifti_image *sinogramImage, nifti_image *bkprImage, nifti_image *attenuationImage, int cam) //FIXME test attenuation
{
    float *sino_data  = (float *) (sinogramImage->data) + cam*bkprImage->nx*bkprImage->ny ;
    float *bkpr_data = (float *)  (bkprImage->data);
    float *attenuation_data = (float *)  (attenuationImage->data);
    float sum_attenuation; 
    int index; 
    for(int y=0; y<bkprImage->ny; y++) {
//#pragma omp parallel for
        for(int x=0; x<bkprImage->nx; x++) {
            sum_attenuation = 0;
            for(int z=0; z<bkprImage->nz; z++) {
                index = x + y*bkprImage->nx + z*bkprImage->nx*bkprImage->ny;
                sum_attenuation += attenuation_data[index]; 
                bkpr_data[index] = sino_data[x + y*bkprImage->nx]*exp(-sum_attenuation);
            }
        }
    }

}


