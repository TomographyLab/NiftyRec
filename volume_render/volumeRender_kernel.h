
//typedef unsigned char VolumeType;
typedef unsigned short VolumeType;
typedef unsigned int  uint;
typedef unsigned char uchar;

extern "C" void setTextureFilterMode(bool bLinearFilter);
extern "C" void initCuda(uint w, uint h, uint d);
extern "C" void freeCudaBuffers();
extern "C" void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint imageW, uint imageH, 
				   float density1, float brightness1, float transferOffset1, float transferScale1, 
                                   float density2, float brightness2, float transferOffset2, float transferScale2);
extern "C" void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix);
extern "C" void setTexture1(VolumeType *h_volume, uint w, uint h, uint d);
extern "C" void setTexture2(VolumeType *h_volume, uint w, uint h, uint d);
extern "C" void setTransferFunction1(float *colormap, unsigned int colormap_size);
extern "C" void setTransferFunction2(float *colormap, unsigned int colormap_size);
extern "C" void setCameraPosition(float x, float y, float z);
extern "C" void setRayStep(float step);
