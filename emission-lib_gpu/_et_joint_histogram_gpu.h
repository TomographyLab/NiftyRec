/*
 *  _et_joint_histogram_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#define BLOCK 256

void et_joint_histogram_gpu(float **d_array_A, float **d_array_B, int **d_jont_hist, int array_size, int hist_size, float min_A, float max_A, float min_B, float max_B);
