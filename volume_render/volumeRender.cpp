/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
 
#include <volumeRender.h>
#include <unistd.h>
#include "_reg_blocksize_gpu.h"

//char *volumeFilename = "test.raw";
uint Vol_width;
uint Vol_height;
uint Vol_depth;
uint width;
uint height;

dim3 blockSize(16, 16);
dim3 gridSize;

float3 viewRotation;
float3 viewTranslation = make_float3(0.0, 0.0, -4.0f);
float invViewMatrix[12];

float density1 = 0.05f;
float brightness1 = 1.0f;
float transferOffset1 = 0.0f;
float transferScale1 = 1.0f;
float density2 = 0.05f;
float brightness2 = 1.0f;
float transferOffset2 = 0.0f;
float transferScale2 = 1.0f;
bool linearFiltering = true;

#define DELTA_STEP 0.0005f
#define DELTA_CAMERA 0.1f
float _step=0.004f;
float _x_camera=0.0f;
float _y_camera=0.0f;
float _z_camera=1000.0f;

GLuint pbo = 0;     // OpenGL pixel buffer object
GLuint tex = 0;     // OpenGL texture object
struct cudaGraphicsResource *cuda_pbo_resource; // CUDA Graphics Resource (to transfer PBO)

int window=0;
//unsigned int timer = 0;
//unsigned int timer_display = 0;

int fpsCount = 0;        // FPS count for averaging
int fpsLimit = 10;       // FPS limit for sampling
int fpsMax   = 20;
bool screenshotRun   = false;
int  screenshotLimit  = 10;
//bool screenshotSingle = false;
//int  screenshotPauseSingle = 100;
unsigned int screenshotCount = 0;
bool has_new_volume1=false;
bool has_new_volume2=false;
bool has_colormap1=false;
bool has_colormap2=false;
bool do_stop=false;
bool stopped=false;
float *h_new_colormap1;
unsigned int h_colormap_size1;
float *h_new_colormap2;
unsigned int h_colormap_size2;
VolumeType *h_new_volume1=NULL;
VolumeType *h_new_volume2=NULL;
VolumeType *image=NULL;
bool do_image_dump=false;

int g_Index = 0;
unsigned int frameCount = 0;
unsigned int g_TotalErrors = 0;
bool g_Verify = false;
bool g_bQAReadback = false;
bool g_bQAGLVerify = false;
bool g_bFBODisplay = false;


void initPixelBuffer();
void stop_now(void);
void close(void);


int loadVolume1()
{
    if (has_new_volume1) {
        setTexture1(h_new_volume1, Vol_width, Vol_height, Vol_depth);
        has_new_volume1=false;
    }
    return 0;
}

int loadVolume2()
{
    if (has_new_volume2) {
        setTexture2(h_new_volume2, Vol_width, Vol_height, Vol_depth);
        has_new_volume2=false;
    }
    return 0;
}

int setColormap1()
{
    if (has_colormap1) {
        setTransferFunction1(h_new_colormap1, h_colormap_size1);
        has_colormap1=false;
    }
    return 0;
}

int setColormap2()
{
    if (has_colormap2) {
        setTransferFunction2(h_new_colormap2, h_colormap_size2);
        has_colormap2=false;
    }
    return 0;
}

int ScreenShot()
{

    screenshotCount++;
    if (screenshotRun==true) {
        if (int(screenshotCount) == screenshotLimit-1) {
            WindowDump(0);
            screenshotCount=0;
        }
    }
    return 0;
}

void computeFPS()
{
    frameCount++;
    fpsCount++;
    if (fpsCount == fpsLimit-1) {
        g_Verify = true;
    }
    if (fpsCount == fpsLimit) {
        //char fps[256];
        //float ifps = 1.f / (gpuGetAverageTimerValue(timer) / 1000.f);
        //sprintf(fps, "Volume Render: %3.1f fps", ifps);  

        //glutSetWindowTitle(fps);
        fpsCount = 0; 
        //CUDA_CHECK_ERROR(cutResetTimer(timer));  
    }
}

// render image using CUDA
void render()
{
    copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);

    // map PBO to get CUDA device pointer
    uint *d_output;
    // map PBO to get CUDA device pointer
    CUDA_SAFE_CALL(cudaGraphicsMapResources(1, &cuda_pbo_resource, 0));
    size_t num_bytes; 
    CUDA_SAFE_CALL(cudaGraphicsResourceGetMappedPointer((void **)&d_output, &num_bytes, cuda_pbo_resource));
    //printf("CUDA mapped PBO: May access %ld bytes\n", num_bytes);

    // clear image
    CUDA_SAFE_CALL(cudaMemset(d_output, 0, width*height*4));

    // call CUDA kernel, writing results to PBO
    render_kernel(gridSize, blockSize, d_output, width, height, density1, brightness1, transferOffset1, transferScale1, 
density2, brightness2, transferOffset2, transferScale2);

    CUDA_CHECK_MSG("kernel failed");

    CUDA_SAFE_CALL(cudaGraphicsUnmapResources(1, &cuda_pbo_resource, 0));
}

// display results using OpenGL (called by GLUT)
void display()
{
    //cutStopTimer(timer_display);
    //float time_since_display = cutGetAverageTimerValue(timer_display);
//    fprintf(stderr,"\ntime: %f",time_since_display);
    //if(not(time_since_display>(1000/fpsMax))) {
//        cutStartTimer(timer_display);
    //    return;
    //}
    //cutResetTimer(timer_display);
    //cutStartTimer(timer_display);

    //CUDA_CHECK_ERROR(cutStartTimer(timer));  

    // use OpenGL to build view matrix
    GLfloat modelView[16];
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotatef(-viewRotation.x, 1.0, 0.0, 0.0);
    glRotatef(-viewRotation.y, 0.0, 1.0, 0.0);
    glTranslatef(-viewTranslation.x, -viewTranslation.y, -viewTranslation.z);
    glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
    glPopMatrix();

    invViewMatrix[0] = modelView[0]; invViewMatrix[1] = modelView[4]; invViewMatrix[2] = modelView[8]; invViewMatrix[3] = modelView[12];
    invViewMatrix[4] = modelView[1]; invViewMatrix[5] = modelView[5]; invViewMatrix[6] = modelView[9]; invViewMatrix[7] = modelView[13];
    invViewMatrix[8] = modelView[2]; invViewMatrix[9] = modelView[6]; invViewMatrix[10] = modelView[10]; invViewMatrix[11] = modelView[14];

    render();

    // display results
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image from PBO
    glDisable(GL_DEPTH_TEST);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

#if 0
    // draw using glDrawPixels (slower)
    glRasterPos2i(0, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
    glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
#else

    // draw using texture

    // copy from pbo to texture
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

    // draw textured quad
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0); glVertex2f(0, 0);
    glTexCoord2f(1, 0); glVertex2f(1, 0);
    glTexCoord2f(1, 1); glVertex2f(1, 1);
    glTexCoord2f(0, 1); glVertex2f(0, 1);
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
#endif
    glutSwapBuffers();
    glutReportErrors();

    //CUDA_CHECK_ERROR(cutStopTimer(timer));  

    computeFPS();
    ScreenShot();

    //cutStopTimer(timer_display);
}

void idle()
{
    if (!stopped)
    {
        glutPostRedisplay();
  
        if (do_stop)
            stop_now();
        loadVolume1();
        loadVolume2();
        setColormap1();
        setColormap2();
  
        if (do_image_dump) {
             WindowDump(1);
             do_image_dump=false;
        }
    }
    usleep(5000);
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key) {
        case 27:
            exit(0);
            break;
        case 'b':
            linearFiltering = !linearFiltering;
			setTextureFilterMode(linearFiltering);
            break;

        case '2':
            density1 += 0.01;
            break;
        case '1':
            density1 -= 0.01;
            break;
        case '"':
            density1 += 0.1;
            break;
        case '!':
            density1 -= 0.1;
            break;
        case 'w':
            brightness1 += 0.1;
            break;
        case 'q':
            brightness1 -= 0.1;
            break;
        case 's':
            transferOffset1 += 0.01;
            break;
        case 'a':
            transferOffset1 -= 0.01;
            break;
        case 'x':
            transferScale1 += 0.01;
            break;
        case 'z':
            transferScale1 -= 0.01;
            break;

        case '4':
            density2 += 0.01;
            break;
        case '3':
            density2 -= 0.01;
            break;
        case 'r':
            brightness2 += 0.1;
            break;
        case 'e':
            brightness2 -= 0.1;
            break;
        case 'f':
            transferOffset2 += 0.01;
            break;
        case 'd':
            transferOffset2 -= 0.01;
            break;
        case 'v':
            transferScale2 += 0.01;
            break;
        case 'c':
            transferScale2 -= 0.01;
            break;

        case 'n':
            screenshotRun=true;
            break;
        case 'm':
            screenshotRun=false;
            break;
        case 't':
            _step += DELTA_STEP;
            setRayStep(_step);
            break;
        case 'y':
            _step -= DELTA_STEP;
            _step=_step<=DELTA_STEP?DELTA_STEP:_step;
            setRayStep(_step);
            break;
        case 'g':
            _z_camera += DELTA_CAMERA;
            setCameraPosition(_x_camera,_y_camera,_z_camera);
            break;
        case 'h':
            _z_camera -= DELTA_CAMERA;
            _z_camera=_z_camera<=DELTA_CAMERA?DELTA_CAMERA:_z_camera;
            setCameraPosition(_x_camera,_y_camera,_z_camera);
            break;


//        case 'w':
//            screeshotSingle=true;
//            break;
        default:
            break;
    }
    fprintf(stderr,"Vol 1\tDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f\nVol 2\tDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f\nGlobal\tRayStep %f zCamera %f\n",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2,_step,_z_camera);
    glutPostRedisplay();
}

int ox, oy;
int buttonState = 0;

void mouse(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
        buttonState  |= 1<<button;
    else if (state == GLUT_UP)
        buttonState = 0;

    ox = x; oy = y;
    glutPostRedisplay();
}

void motion(int x, int y)
{
    float dx, dy;
    dx = x - ox;
    dy = y - oy;

    if (buttonState == 4) {
        // right = zoom
        viewTranslation.z += 100*(dy / 100.0);
    } 
    else if (buttonState == 2) {
        // middle = translate
        viewTranslation.x += dx / 100.0;
        viewTranslation.y -= dy / 100.0;
    }
    else if (buttonState == 1) {
        // left = rotate
        viewRotation.x += dy / 5.0;
        viewRotation.y += dx / 5.0;
    }

    ox = x; oy = y;
    glutPostRedisplay();
}

int iDivUp(int a, int b){
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

void reshape(int w, int h)
{
    width = w; height = h;
    initPixelBuffer();

    // calculate new grid size
    gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y));

    glViewport(0, 0, w, h);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
}

void cleanup()
{
    //CUDA_CHECK_ERROR( cutDeleteTimer( timer));
    freeCudaBuffers();
    if (pbo) {
	//CUDA_SAFE_CALL(cudaGraphicsUnregisterResource(cuda_pbo_resource));
        glDeleteBuffers(1, &pbo);
        glDeleteTextures(1, &tex);
    }
    free(h_new_volume1);
    free(h_new_volume2);
    free(h_new_colormap1);
    free(h_new_colormap2);
    free(image);
    cudaThreadExit();
}

int initGL(int *argc, char **argv)
{
    // initialize GLUT callback functions
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    window = glutCreateWindow("CUDA volume rendering");

    glewInit();
    if (!glewIsSupported("GL_VERSION_2_0")) {
        fprintf(stderr,"\nRequired OpenGL extensions missing.");
        return 1;
    }
    return 0;
}

void initPixelBuffer()
{
/*    if (pbo) {
        // unregister this buffer object from CUDA C
        CUDA_SAFE_CALL(cudaGraphicsUnregisterResource(cuda_pbo_resource));

        // delete old buffer
        glDeleteBuffers(1, &pbo);
        glDeleteTextures(1, &tex);
    }*/

    // create pixel buffer object for display
    glGenBuffers(1, &pbo);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
    glBufferData(GL_PIXEL_UNPACK_BUFFER, width*height*sizeof(GLubyte)*4, 0, GL_STREAM_DRAW);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

    // register this buffer object with CUDA
    CUDA_SAFE_CALL(cudaGraphicsGLRegisterBuffer(&cuda_pbo_resource, pbo, cudaGraphicsMapFlagsWriteDiscard));	

    // create texture for display
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Load raw data from disk
void *loadRawFile(char *filename, size_t size)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }

	void *data = malloc(size);
//	size_t read = fread(data, 1, size, fp);
	fclose(fp);

    return data;
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
extern "C" void run_gui(int volume_x, int volume_y, int volume_z, int camera_x, int camera_y) 
{
    // First initialize OpenGL context, so we can properly set the GL for CUDA.
    // This is necessary in order to achieve optimal performance with OpenGL/CUDA interop.

    Vol_width  = volume_x;
    Vol_height = volume_y;
    Vol_depth  = volume_z;
    width = camera_x;
    height = camera_y;

    int argc=0;
    char **argv=NULL;

    cudaGLSetGLDevice( gpuGetMaxGflopsDeviceId() );
    if (initGL( &argc, argv )) return;

    // load volume data
    size_t size = Vol_width*Vol_height*Vol_depth*sizeof(VolumeType);
    void *h_volume1 = malloc(size); //loadRawFile(volumeFilename, size);
    void *h_volume2 = malloc(size);
    memset(h_volume2,0,size);

    initCuda(Vol_width, Vol_height, Vol_depth); 
    setTexture1((VolumeType*)h_volume1,Vol_width, Vol_height, Vol_depth);
    setTexture2((VolumeType*)h_volume2,Vol_width, Vol_height, Vol_depth);
    free(h_volume1);
    free(h_volume2);
    h_new_volume1 = (VolumeType*)malloc(Vol_width*Vol_height*Vol_depth*sizeof(VolumeType));
    h_new_volume2 = (VolumeType*)malloc(Vol_width*Vol_height*Vol_depth*sizeof(VolumeType));
   
    setCameraPosition(_x_camera, _y_camera, _z_camera);
    setRayStep(_step);

    //CUDA_CHECK_ERROR( cutCreateTimer( &timer));
    //CUDA_CHECK_ERROR( cutCreateTimer( &timer_display));
    //cutStartTimer(timer_display);
 
    // calculate new grid size
    gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y));

    // This is the normal rendering path for VolumeRender
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    //glutCloseFunc(close);
    initPixelBuffer();
   // glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
   // glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
//    atexit(cleanup);

    do_stop=false;
    stopped=false;
    glutMainLoop();
}


int WindowDump(bool single)
{
   int i,j;
   FILE *fptr;
   static int counter = 0; /* This supports animation sequences */
   char fname[256];
   unsigned char *image;
   int width =  glutGet(GLUT_WINDOW_WIDTH);
   int height = glutGet(GLUT_WINDOW_WIDTH);

   /* Allocate our buffer for the image */
   if ((image = (unsigned char*)malloc(3*width*height*sizeof(char))) == NULL) {
      fprintf(stderr,"Failed to allocate memory for image\n");
      return(false);
   }
   glPixelStorei(GL_PACK_ALIGNMENT,1);

   /* Open the file */

   if (single)
       sprintf(fname,"./screenshot.raw");
   else
       sprintf(fname,"./C_%04d.raw",counter);
   if ((fptr = fopen(fname,"w")) == NULL) {
      fprintf(stderr,"Failed to open file for window dump\n");
      return(false);
   }

   /* Copy the image into our buffer */
   glReadBuffer(GL_BACK_LEFT);
   glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);

   /* Write the raw file */
   /* fprintf(fptr,"P6\n%d %d\n255\n",width,height); for ppm */
   for (j=height-1;j>=0;j--) {
      for (i=0;i<width;i++) {
         fputc(image[3*j*width+3*i+0],fptr);
         fputc(image[3*j*width+3*i+1],fptr);
         fputc(image[3*j*width+3*i+2],fptr);
      }
   }
   fclose(fptr);

   if (!single)
       counter++;

   /* Clean up */
   free(image);
   return(true);

}




extern "C" void stop(void)
{
    do_stop=true;
}

void stop_now(void)
{
    stopped=true;
    do_stop=false;
    cleanup();
    glutDestroyWindow(window);
}

void close(void)
{
//    stopped=true;
//    fprintf(stderr,"\nclose");
//    cleanup();
//    glutDestroyWindow(window);
//    while(1);
}

extern "C" int get_density1(void)
{
    return density1;
}

extern "C" int set_density1(float _density)
{
    density1 = _density;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return density1;
}

extern "C" int get_brightness1(void)
{
    return brightness1;
}

extern "C" int set_brightness1(float _brightness)
{
    brightness1 = _brightness;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return brightness1;
}

extern "C" int get_transferOffset1(void)
{
    return transferOffset1;
}

extern "C" int set_transferOffset1(float _transferOffset)
{
    transferOffset1 = _transferOffset;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return transferOffset1;
}

extern "C" int get_transferScale1(void)
{
    return transferScale1;
}

extern "C" int set_transferScale1(float _transferScale)
{
    transferScale1 = _transferScale;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return transferScale1;
}

extern "C" int set_volume1(VolumeType *volume)
{
    memcpy(h_new_volume1,volume,Vol_width*Vol_height*Vol_depth*sizeof(VolumeType));
    has_new_volume1=true;
    return 0;
}

extern "C" int get_density2(void)
{
    return density2;
}

extern "C" int set_density2(float _density)
{
    density2 = _density;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return density2;
}

extern "C" int get_brightness2(void)
{
    return brightness2;
}

extern "C" int set_brightness2(float _brightness)
{
    brightness2 = _brightness;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return brightness2;
}

extern "C" int get_transferOffset2(void)
{
    return transferOffset2;
}

extern "C" int set_transferOffset2(float _transferOffset)
{
    transferOffset2 = _transferOffset;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return transferOffset2;
}

extern "C" int get_transferScale2(void)
{
    return transferScale2;
}

extern "C" int set_transferScale2(float _transferScale)
{
    transferScale2 = _transferScale;
    glutPostRedisplay();
    fprintf(stderr,"\nDens1: %f  Brigh1: %f  Offs1: %f  Scale1: %f     \nDens2: %f  Brigh2: %f  Offs2: %f  Scale2: %f",
            density1,brightness1,transferOffset1,transferScale1,density2,brightness2,transferOffset2,transferScale2);
    return transferScale2;
}

extern "C" int set_volume2(VolumeType *volume)
{
    memcpy(h_new_volume2,volume,Vol_width*Vol_height*Vol_depth*sizeof(VolumeType));
    has_new_volume2=true;
    return 0;
}

extern "C" int set_colormap1(float *colormap, unsigned int colormap_size)
{
    h_new_colormap1 = (float*)realloc((void*)h_new_colormap1,colormap_size*4*sizeof(float));
    memcpy(h_new_colormap1,colormap,colormap_size*4*sizeof(float));
    has_colormap1=true;
    h_colormap_size1 = colormap_size;
    return 0;
}

extern "C" int set_colormap2(float *colormap, unsigned int colormap_size)
{
    h_new_colormap2 = (float*)realloc((void*)h_new_colormap2,colormap_size*4*sizeof(float));
    memcpy(h_new_colormap2,colormap,colormap_size*4*sizeof(float));
    h_colormap_size2 = colormap_size;
    has_colormap2=true;
    return 0;
}



/*extern "C" unsigned char *get_screenshot(uint *size_x, uint *size_y)
{
   int width =  glutGet(GLUT_WINDOW_WIDTH);
   int height = glutGet(GLUT_WINDOW_WIDTH);
   if ((image = (unsigned char*)realloc(image, 3*width*height*sizeof(unsigned char))) == NULL) {
      fprintf(stderr,"Failed to allocate memory for image\n");
      return(false);
   }
   glPixelStorei(GL_PACK_ALIGNMENT,1);
   glReadBuffer(GL_BACK_LEFT);
   glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
   memset(image,100,3*width*height*sizeof(unsigned char));

   *size_x = width;
   *size_y = height;
   return image;
}*/

extern "C" void dump_screenshot()
{
    do_image_dump=true;
    while(do_image_dump)
    {
        usleep(1000);
    }
}


extern "C" void translate(float x, float y)
{
    viewTranslation.x += x;
    viewTranslation.y += y;
    glutPostRedisplay();
}

extern "C" void rotate(float x, float y)
{
    viewRotation.x += x;
    viewRotation.y += y;
    glutPostRedisplay();
}

extern "C" void zoom(float z)
{
    viewTranslation.z += z;
    glutPostRedisplay();
}
