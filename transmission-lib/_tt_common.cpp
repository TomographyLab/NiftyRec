/*
 *  _tt_common.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include <_tt_common.h>
#include <stdlib.h>

mat_44 reg_mat_44_mul(mat_44 *A, mat_44 *B)
{
	mat_44 R;
	
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			R.m[i][j] = A->m[i][0]*B->m[0][j] + A->m[i][1]*B->m[1][j] + A->m[i][2]*B->m[2][j] + A->m[i][3]*B->m[3][j];
		}
	}
	
	return R;
}


int create_rotation_matrix44(mat_44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z)
{
	int status = 0;
	float s_theta_x, c_theta_x, s_theta_y, c_theta_y, s_theta_z, c_theta_z;
	mat_44 *rotation_x = (mat_44 *)calloc(1,sizeof(mat_44));
	mat_44 *rotation_y = (mat_44 *)calloc(1,sizeof(mat_44));
	mat_44 *rotation_z = (mat_44 *)calloc(1,sizeof(mat_44));
	
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
	rotation_z->m[0][3] = - center_x*c_theta_z + center_y*s_theta_z + center_x;
	rotation_z->m[1][0] = s_theta_z;
	rotation_z->m[1][1] = c_theta_z;
	rotation_z->m[1][2] = 0.0;
	rotation_z->m[1][3] = - center_x*s_theta_z - center_y*c_theta_z + center_y;
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
	rotation_y->m[0][3] = - c_theta_y*center_x - s_theta_y*center_z + center_x;
	rotation_y->m[1][0] = 0.0;
	rotation_y->m[1][1] = 1.0;
	rotation_y->m[1][2] = 0.0;
	rotation_y->m[1][3] = 0.0;
	rotation_y->m[2][0] = -s_theta_y;
	rotation_y->m[2][1] = 0.0;
	rotation_y->m[2][2] = c_theta_y;
	rotation_y->m[2][3] = s_theta_y*center_x - c_theta_y*center_z + center_z;
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
	rotation_x->m[1][3] = -c_theta_x*center_y + s_theta_x*center_z + center_y;
	rotation_x->m[2][0] = 0.0;
	rotation_x->m[2][1] = s_theta_x;
	rotation_x->m[2][2] = c_theta_x;
	rotation_x->m[2][3] = -s_theta_x*center_y - c_theta_x*center_z + center_z;
	rotation_x->m[3][0] = 0.0;
	rotation_x->m[3][1] = 0.0;
	rotation_x->m[3][2] = 0.0;
	rotation_x->m[3][3] = 1.0;	
	
	*transformationMatrix = reg_mat_44_mul(rotation_y, rotation_x);
	*transformationMatrix = reg_mat_44_mul(rotation_z, transformationMatrix);
	
	free(rotation_x);
	free(rotation_y);
	free(rotation_z);

	return status;
}

