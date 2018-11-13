
#include <GL/glew.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

// CUDA Runtime, Interop, and includes
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>

// CUDA utilities
#include <helper_cuda.h>
#include <helper_cuda_gl.h>

// Helper functions
#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_timer.h>

#include <volumeRender_kernel.h>

typedef unsigned int uint;
typedef unsigned char uchar;

extern "C" void run_gui(int volume_x, int volume_y, int volume_z, int camera_x, int camera_y);

extern "C" int WindowDump(bool single);
extern "C" int ScreenShot(void);


#define MAX_EPSILON_ERROR 5.00f
#define THRESHOLD         0.30f
