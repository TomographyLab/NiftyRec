/*
 *  _pet_line_integrals_compressed_gpu.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Harvard University, Martinos Center for Biomedical Imaging
 *  Jan. 2014. 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "_reg_blocksize_gpu.h"
#include "nifti1_io.h"

// note: block_size affects the performance, a good value is 512 
    
void pet_line_integral_compressed_gpu(float *d_activity, float *d_attenuation, float *d_projection, unsigned short *d_locations, unsigned int N_locations, unsigned int N_u, unsigned int N_v, unsigned int N_samples, unsigned int direction, unsigned int block_size); 
