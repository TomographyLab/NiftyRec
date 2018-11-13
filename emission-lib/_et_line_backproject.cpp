/*
 *  _et_line_backproject.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et_line_backproject.h"
#ifdef _OPENMP
#include "omp.h"
#endif

void et_line_backproject(nifti_image *sinogramImage, nifti_image *bkprImage, int cam)
{
    float *sino_data  = (float *) (sinogramImage->data) + cam*bkprImage->nx*bkprImage->ny ;
    float *bkpr_data = (float *)  (bkprImage->data);
    for(int y=0; y<bkprImage->ny; y++) {
//#pragma omp parallel for
        for(int x=0; x<bkprImage->nx; x++) {
            for(int z=0; z<bkprImage->nz; z++) {
                bkpr_data[x + y*bkprImage->nx + z*bkprImage->nx*bkprImage->ny] = sino_data[x + y*bkprImage->nx];
            }
        }
    }

}


