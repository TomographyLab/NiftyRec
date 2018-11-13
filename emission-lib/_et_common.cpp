/*
 *  _et_common.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et_common.h"


extern "C" int fprintf_verbose (const char *__restrict __format, ...)
{
    va_list args;
    va_start(args,__format);  
    int return_status = 1;

    #ifdef _VERBOSE
    return_status = vfprintf(stderr, __format, args);
    #else
    return_status = 0;
    #endif

     va_end(args);
     return return_status;
 }





int et_create_rotation_matrix(mat44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, int axis_order)
{
	int status = 0;
	float s_theta_x, c_theta_x, s_theta_y, c_theta_y, s_theta_z, c_theta_z;
	mat44 *rotation_x = (mat44 *)calloc(1,sizeof(mat44));
	mat44 *rotation_y = (mat44 *)calloc(1,sizeof(mat44));
	mat44 *rotation_z = (mat44 *)calloc(1,sizeof(mat44));
	mat44 *translation = (mat44 *)calloc(1,sizeof(mat44));
	
	// Initialize affine transform matrix
	s_theta_x = sin(theta_x);
	c_theta_x = cos(theta_x);
	s_theta_y  = sin(theta_y);
	c_theta_y  = cos(theta_y);
	s_theta_z = sin(theta_z);
	c_theta_z = cos(theta_z);
	
	rotation_z->m[0][0] = c_theta_z;
	rotation_z->m[0][1] = -s_theta_z;
	rotation_z->m[0][2] = 0.0;
	rotation_z->m[0][3] = 0.0; 
	rotation_z->m[1][0] = s_theta_z;
	rotation_z->m[1][1] = c_theta_z;
	rotation_z->m[1][2] = 0.0;
	rotation_z->m[1][3] = 0.0;
	rotation_z->m[2][0] = 0.0;
	rotation_z->m[2][1] = 0.0;
	rotation_z->m[2][2] = 1.0;
	rotation_z->m[2][3] = 0.0;
	rotation_z->m[3][0] = 0.0;
	rotation_z->m[3][1] = 0.0;
	rotation_z->m[3][2] = 0.0;
	rotation_z->m[3][3] = 1.0;	
	
	rotation_y->m[0][0] = c_theta_y;
	rotation_y->m[0][1] = 0.0;
	rotation_y->m[0][2] = s_theta_y;
	rotation_y->m[0][3] = 0.0;
	rotation_y->m[1][0] = 0.0;
	rotation_y->m[1][1] = 1.0;
	rotation_y->m[1][2] = 0.0;
	rotation_y->m[1][3] = 0.0;
	rotation_y->m[2][0] = -s_theta_y;
	rotation_y->m[2][1] = 0.0;
	rotation_y->m[2][2] = c_theta_y;
	rotation_y->m[2][3] = 0.0; 
	rotation_y->m[3][0] = 0.0;
	rotation_y->m[3][1] = 0.0;
	rotation_y->m[3][2] = 0.0;
	rotation_y->m[3][3] = 1.0;	
	
	rotation_x->m[0][0] = 1.0;
	rotation_x->m[0][1] = 0.0;
	rotation_x->m[0][2] = 0.0;
	rotation_x->m[0][3] = 0.0;
	rotation_x->m[1][0] = 0.0;
	rotation_x->m[1][1] = c_theta_x;
	rotation_x->m[1][2] = -s_theta_x;
	rotation_x->m[1][3] = 0.0; 
	rotation_x->m[2][0] = 0.0;
	rotation_x->m[2][1] = s_theta_x;
	rotation_x->m[2][2] = c_theta_x;
	rotation_x->m[2][3] = 0.0; 
	rotation_x->m[3][0] = 0.0;
	rotation_x->m[3][1] = 0.0;
	rotation_x->m[3][2] = 0.0;
	rotation_x->m[3][3] = 1.0;	
	
	translation->m[0][0] = 0.0; 
	translation->m[0][1] = 0.0;
	translation->m[0][2] = 0.0;
	translation->m[0][3] = center_x;
	translation->m[1][0] = 0.0; 
	translation->m[1][1] = 0.0;
	translation->m[1][2] = 0.0;
	translation->m[1][3] = center_y;
	translation->m[2][0] = 0.0; 
	translation->m[2][1] = 0.0;
	translation->m[2][2] = 0.0;
	translation->m[2][3] = center_z;
	translation->m[3][0] = 0.0; 
	translation->m[3][1] = 0.0;
	translation->m[3][2] = 0.0;
	translation->m[3][3] = 1;
	
	switch (axis_order) {
		case XYZ_ROTATION:
			*transformationMatrix = reg_mat44_mul(rotation_y, rotation_x);
			*transformationMatrix = reg_mat44_mul(rotation_z, transformationMatrix);			
			break;
		case XZY_ROTATION:
			*transformationMatrix = reg_mat44_mul(rotation_z, rotation_x);
			*transformationMatrix = reg_mat44_mul(rotation_y, transformationMatrix);			
			break;
		case YXZ_ROTATION:
			*transformationMatrix = reg_mat44_mul(rotation_x, rotation_y);
			*transformationMatrix = reg_mat44_mul(rotation_z, transformationMatrix);			
			break;
		case YZX_ROTATION:
			*transformationMatrix = reg_mat44_mul(rotation_z, rotation_y);
			*transformationMatrix = reg_mat44_mul(rotation_x, transformationMatrix);			
			break;
		case ZXY_ROTATION:
			*transformationMatrix = reg_mat44_mul(rotation_x, rotation_z);
			*transformationMatrix = reg_mat44_mul(rotation_y, transformationMatrix);			
			break;
		case ZYX_ROTATION:
			*transformationMatrix = reg_mat44_mul(rotation_y, rotation_z);
			*transformationMatrix = reg_mat44_mul(rotation_x, transformationMatrix);			
			break;
		default: //default: XYZ_ROTATION
			*transformationMatrix = reg_mat44_mul(rotation_y, rotation_x);
			*transformationMatrix = reg_mat44_mul(rotation_z, transformationMatrix);
			break;
	}
	
	
	*translation = reg_mat44_mul(transformationMatrix, translation);
	
	transformationMatrix->m[0][3] = center_x - translation->m[0][3];
	transformationMatrix->m[1][3] = center_y - translation->m[1][3];
	transformationMatrix->m[2][3] = center_z - translation->m[2][3];
	
	free(rotation_x);
	free(rotation_y);
	free(rotation_z);
    free(translation); 

//	fprintf(stderr,"=============== Input: ==============\n");	
//	fprintf(stderr,"Cx: %3.3f  Cy: %3.3f  Cz: %3.3f  \n",center_x, center_y, center_z); 
//	fprintf(stderr,"Rx: %3.3f  Ry: %3.3f  Rz: %3.3f  \n",theta_x, theta_y, theta_z); 	
//	fprintf(stderr,"=======    Rotation matrix:    ======\n");
//	for (int i=0; i<4; i++)
//			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",transformationMatrix->m[i][0],transformationMatrix->m[i][1],transformationMatrix->m[i][2],transformationMatrix->m[i][3]);
//    fprintf(stderr,"=====================================\n"); 
	return status;
}


int et_create_scale_matrix(mat44 *transformationMatrix, float scale_x, float scale_y, float scale_z)
{
	int status = 0;
	
	transformationMatrix->m[0][0] = scale_x;
	transformationMatrix->m[0][1] = 0.0;
	transformationMatrix->m[0][2] = 0.0;
	transformationMatrix->m[0][3] = 0.0; 
	transformationMatrix->m[1][0] = 0.0;
	transformationMatrix->m[1][1] = scale_y;
	transformationMatrix->m[1][2] = 0.0;
	transformationMatrix->m[1][3] = 0.0;
	transformationMatrix->m[2][0] = 0.0;
	transformationMatrix->m[2][1] = 0.0;
	transformationMatrix->m[2][2] = scale_z;
	transformationMatrix->m[2][3] = 0.0;
	transformationMatrix->m[3][0] = 0.0;
	transformationMatrix->m[3][1] = 0.0;
	transformationMatrix->m[3][2] = 0.0;
	transformationMatrix->m[3][3] = 1.0;	

/*
	fprintf(stderr,"=============== Input: ==============\n");	
	fprintf(stderr,"Sx: %3.3f  Sy: %3.3f  Sz: %3.3f  \n",scale_x, scale_y, scale_z);  	
	fprintf(stderr,"=======     Scale matrix:      ======\n");
	for (int i=0; i<4; i++)
			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",transformationMatrix->m[i][0],transformationMatrix->m[i][1],transformationMatrix->m[i][2],transformationMatrix->m[i][3]);
    fprintf(stderr,"=====================================\n"); 
*/

	return status;
}


int et_create_translation_matrix(mat44 *transformationMatrix, float t_x, float t_y, float t_z)
{
	int status = 0;
	
	transformationMatrix->m[0][0] = 1.0;
	transformationMatrix->m[0][1] = 0.0;
	transformationMatrix->m[0][2] = 0.0;
	transformationMatrix->m[0][3] = t_x; 
	transformationMatrix->m[1][0] = 0.0;
	transformationMatrix->m[1][1] = 1.0;
	transformationMatrix->m[1][2] = 0.0;
	transformationMatrix->m[1][3] = t_y;
	transformationMatrix->m[2][0] = 0.0;
	transformationMatrix->m[2][1] = 0.0;
	transformationMatrix->m[2][2] = 1.0;
	transformationMatrix->m[2][3] = t_z;
	transformationMatrix->m[3][0] = 0.0;
	transformationMatrix->m[3][1] = 0.0;
	transformationMatrix->m[3][2] = 0.0;
	transformationMatrix->m[3][3] = 1.0;	

/*
	fprintf(stderr,"=============== Input: ==============\n");	
	fprintf(stderr,"Tx: %3.3f  Ty: %3.3f  Tz: %3.3f  \n",t_x,t_y,t_z);  	
	fprintf(stderr,"=======  Translation matrix:   ======\n");
	for (int i=0; i<4; i++)
			fprintf(stderr,"[%3.3f  %3.3f  %3.3f  %3.3f]\n",transformationMatrix->m[i][0],transformationMatrix->m[i][1],transformationMatrix->m[i][2],transformationMatrix->m[i][3]);
    fprintf(stderr,"=====================================\n"); 
*/

	return status;
}


int et_create_affine_matrix_from_buffer(mat44 *transformationMatrix, float *buffer)
{
	int status = 0;
	
	transformationMatrix->m[0][0] = buffer[0];
	transformationMatrix->m[0][1] = buffer[1];
	transformationMatrix->m[0][2] = buffer[2];
	transformationMatrix->m[0][3] = buffer[3]; 
	transformationMatrix->m[1][0] = buffer[4];
	transformationMatrix->m[1][1] = buffer[5];
	transformationMatrix->m[1][2] = buffer[6];
	transformationMatrix->m[1][3] = buffer[7];
	transformationMatrix->m[2][0] = buffer[8];
	transformationMatrix->m[2][1] = buffer[9];
	transformationMatrix->m[2][2] = buffer[10];
	transformationMatrix->m[2][3] = buffer[11];
	transformationMatrix->m[3][0] = buffer[12];
	transformationMatrix->m[3][1] = buffer[13];
	transformationMatrix->m[3][2] = buffer[14];
	transformationMatrix->m[3][3] = buffer[15];	

	return status;
}




extern "C" char *niftyrec_error_msg(int status)
{

    switch (status)
        {
        case niftyrec_success:
            sprintf(_niftyrec_msg,"niftyrec: OK.\n");
            break;
        case niftyrec_error_unspecified:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - unspecified error.\n");
            break;
        case niftyrec_error_parameters:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - inconsisten paramenters.\n");
            break;
        case niftyrec_error_kernel:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - gpu kernel execution error.\n");
            break;
        case niftyrec_error_allocgpu:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - GPU memory allocation.\n");
            break;
        case niftyrec_error_alloccpu:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - cpu memory allocation.\n");
            break;
        case niftyrec_error_transfergpu:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - transfering memory from guest to GPU or vice versa.\n");
            break;
        case niftyrec_error_nogpubuilt:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - niftyrec was built without GPU support. Please rebuild or install binary built with GPU support in order to use the GPU.\n");
            break;
        case niftyrec_error_nogpuinstalled:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - no GPUs installed in the system.\n");
            break;
        default: 
            sprintf(_niftyrec_msg,"niftyrec: unknown status.\n");
            break;
        }
    return _niftyrec_msg;
}


