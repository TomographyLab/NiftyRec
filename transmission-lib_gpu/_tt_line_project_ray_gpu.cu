/*
 *  _tt_line_project_ray_gpu.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _TTPROJECTRAY_CU_
#define _TTPROJECTRAY_CU_

// Utilities and System includes

//#include <_tt_project_ray_gpu.h>

#define MAX_EPSILON_ERROR 5.00f
#define THRESHOLD         0.30f

#define MAX(a,b) ((a > b) ? a : b)


//#include <cutil_inline.h>
#include "_reg_blocksize_gpu.h"
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
//#include <sys/time.h>
#include <_tt_common.h>


inline int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

int set_inViewMatrix(float *invViewMatrix, float_2 detector_scale, float_3 detector_transl, float_3 detector_rotat)
{
    memset((void*)invViewMatrix,0,12*sizeof(float));
    //rotate
    mat_44 *rotation = (mat_44 *)calloc(1,sizeof(mat_44));
    create_rotation_matrix44(rotation, detector_rotat.x,detector_rotat.y,detector_rotat.z,0,0,0);
    //scale
    mat_44 *scale = (mat_44 *)calloc(1,sizeof(mat_44));
    scale->m[0][0] =detector_scale.x;
    scale->m[1][1] =detector_scale.y;
    scale->m[2][2] =1;
    //transform
    mat_44 *m = (mat_44 *)calloc(1,sizeof(mat_44));
    *m = reg_mat_44_mul(rotation,scale);
    invViewMatrix[0]=m->m[0][0]; invViewMatrix[1]=m->m[0][1]; invViewMatrix[2] =m->m[0][2]; 
    invViewMatrix[4]=m->m[1][0]; invViewMatrix[5]=m->m[1][1]; invViewMatrix[6] =m->m[1][2]; 
    invViewMatrix[8]=m->m[2][0]; invViewMatrix[9]=m->m[2][1]; invViewMatrix[10]=m->m[2][2];
    //translate
    invViewMatrix[3] =detector_transl.x;
    invViewMatrix[7] =detector_transl.y;
    invViewMatrix[11]=detector_transl.z; 
    //cleanup
    free(rotation);
    free(scale);
    free(m);
    return 0;
}


#endif

