
// Utilities and System includes
#include <cutil_inline.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <tt_project_ray_gpu.h>
#include <_tt_common.h>
#include <sys/time.h>

#define MAX_EPSILON_ERROR 5.00f
#define THRESHOLD         0.30f

#define MAX(a,b) ((a > b) ? a : b)

extern "C" void setTextureFilterMode(bool bLinearFilter);
extern "C" void initCuda(void *h_volume, cudaExtent volumeSize);
extern "C" void freeCudaBuffers();
extern "C" void tt_project_ray_kernel(dim3 gridSize, dim3 blockSize, float *d_output, float3 source_position, float3 volume_size, u_int imageW, u_int imageH, float t_step);
extern "C" void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix);


int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

int set_inViewMatrix(float *invViewMatrix, float_2 detector_scale, float_3 detector_transl, float_3 detector_rotat)
{
    memset((void*)invViewMatrix,0,12*sizeof(float));
    //rotate
    mat44 *rotation = (mat44 *)calloc(1,sizeof(mat44));
    create_rotation_matrix44(rotation, detector_rotat.x,detector_rotat.y,detector_rotat.z,0,0,0);
    //scale
    mat44 *scale = (mat44 *)calloc(1,sizeof(mat44));
    scale->m[0][0] =detector_scale.x;
    scale->m[1][1] =detector_scale.y;
    scale->m[2][2] =1;
    //transform
    mat44 *m = (mat44 *)calloc(1,sizeof(mat44));
    *m = reg_mat44_mul(rotation,scale);
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


////////////////////////////////////////////////////////////////////////////////
// tt_project_ray_array
////////////////////////////////////////////////////////////////////////////////
extern "C" int tt_project_ray_array(VolumeType h_volume[], u_int_3 volume_voxels, float out_projections[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_position[], float_3 volume_size, float t_step)
{
    dim3 blockSize(16, 16);
    bool linearFiltering = true;
    dim3 gridSize;

    float  invViewMatrix[12];
    float3 sourcePosition;
    float3 volumeSize;
    float *d_output;
    float *d_output_proj;

    cuInit(0);

    cudaExtent vsize = make_cudaExtent(volume_voxels.x,volume_voxels.y,volume_voxels.z);
fprintf(stderr, "\n %d, %d, %d",volume_voxels.x,volume_voxels.y,volume_voxels.z);
    initCuda(h_volume, vsize);
    setTextureFilterMode(linearFiltering);
    gridSize = dim3(iDivUp(detector_pixels.w, blockSize.x), iDivUp(detector_pixels.h, blockSize.y));

    //Allocate memory for projections on the device
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_output, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float) ));
    CUDA_SAFE_CALL(cudaMemset((void *)d_output,0, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float) ));

    struct timeval start_time; gettimeofday( &start_time, 0);
    struct timeval t_time;
    float elapsed_time;

    volumeSize.x = volume_size.x;
    volumeSize.y = volume_size.y;
    volumeSize.z = volume_size.z;

    for (int proj=0;proj<n_projections;proj++)
    {
        //define invViewMatrix (position of detector) and position of source
        set_inViewMatrix(invViewMatrix, detector_scale[proj], detector_transl[proj], detector_rotat[proj]);
        sourcePosition.x = source_position[proj].x;
        sourcePosition.y = source_position[proj].y;
        sourcePosition.z = source_position[proj].z;

        //project
        copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);
        d_output_proj = (float*) d_output + proj * detector_pixels.w * detector_pixels.h;
        tt_project_ray_kernel(gridSize, blockSize, d_output_proj, sourcePosition, volumeSize, detector_pixels.w, detector_pixels.h, t_step);
    }

    gettimeofday( &t_time, 0);
    elapsed_time = (float) (1000.0 * ( t_time.tv_sec - start_time.tv_sec) + (0.001 * (t_time.tv_usec - start_time.tv_usec)) );
    fprintf(stderr,"\nTime per projection %d %d %d -> %d %d: %f ms",volume_voxels.x,volume_voxels.y,volume_voxels.z,detector_pixels.w,detector_pixels.h,elapsed_time/n_projections);

    //Copy result back to host
    CUDA_SAFE_CALL(cudaMemcpy(out_projections, d_output, n_projections*detector_pixels.w*detector_pixels.h*sizeof(float), cudaMemcpyDeviceToHost));

    //Clean up
    CUDA_SAFE_CALL(cudaFree(d_output));
    freeCudaBuffers();
    cudaThreadExit();
    return 0;
}

