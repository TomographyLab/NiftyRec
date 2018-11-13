
#include <volumeRender.h>
#define DEF_VOLUME_SIZE 512
#define DEF_CAMERA_SIZE 512

int main(int argc, char **argv)
{
    uint volume_size=DEF_VOLUME_SIZE;
    uint camera_size=DEF_CAMERA_SIZE;
    run_gui(volume_size,volume_size,volume_size,camera_size,camera_size);
}
