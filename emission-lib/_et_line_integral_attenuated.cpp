/*
 *  _et_line_integral_attenuated.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et_line_integral_attenuated.h"

void et_line_integral_attenuated(nifti_image *activityImage, nifti_image *attenuationImage, nifti_image *sinoImage, int cam, float background_activity) 
{
    float *sino_data = (float *) (sinoImage->data) + cam*activityImage->nx*activityImage->ny ; 
    float *activity_data = NULL; 
    float *attenuation_data = NULL; 

    float sum_attenuation;
    int index;
    memset((void*) sino_data,0,activityImage->nx*activityImage->ny*sizeof(float));

    if (activityImage!=NULL && attenuationImage!=NULL)
        {
        activity_data = (float *) (activityImage->data); 
        attenuation_data = (float *) (attenuationImage->data); 
        for(int y=0; y<activityImage->ny; y++) {
            for(int x=0; x<activityImage->nx; x++) {
                sum_attenuation = 0;
                for(int z=0; z<activityImage->nz; z++) {
                    index = x + y*activityImage->nx + z*activityImage->nx*activityImage->ny; 
                    sum_attenuation += attenuation_data[index]; 
                    sino_data[x + y*activityImage->nx] = sino_data[x + y*activityImage->nx] + activity_data[index]*exp(-sum_attenuation); 
                    }
                sino_data[x + y*attenuationImage->nx] = sino_data[x + y*attenuationImage->nx] + background_activity * exp(-sum_attenuation); 
                }
            }
        }
    else if (attenuationImage==NULL)
        {
        activity_data = (float *) (activityImage->data); 
        for(int y=0; y<activityImage->ny; y++) {
            for(int x=0; x<activityImage->nx; x++) {
                for(int z=0; z<activityImage->nz; z++) {
                    index = x + y*activityImage->nx + z*activityImage->nx*activityImage->ny;  
                    sino_data[x + y*activityImage->nx] = sino_data[x + y*activityImage->nx] + activity_data[index]; 
                    }
                sino_data[x + y*activityImage->nx] = sino_data[x + y*activityImage->nx] + background_activity;  
                }
            }
        }
    else 
        {
        attenuation_data = (float *) (attenuationImage->data); 
        for(int y=0; y<attenuationImage->ny; y++) {
            for(int x=0; x<attenuationImage->nx; x++) {
                sum_attenuation = 0;
                for(int z=0; z<attenuationImage->nz; z++) {
                    index = x + y*attenuationImage->nx + z*attenuationImage->nx*attenuationImage->ny; 
                    sum_attenuation += attenuation_data[index]; 
                    }
                sino_data[x + y*attenuationImage->nx] = background_activity * exp(-sum_attenuation); 
                }
            }
        }

}


