
from volumeRender import VolumeRender
from numpy import ones, zeros, uint16, double
from time import sleep
import os

N = 128
thickness = 2
image_steps = 15
dynrange = 4.5

def run():
    if not 'screenshots' in os.listdir('./'):
        os.mkdir('./screenshots')
    #load Matlab data
#    import scipy.io
#    mri = scipy.io.loadmat('brainweb_128.mat')
#    activity = scipy.io.loadmat('L.mat')
#    v1 = uint16( activity['L']*(2**16*(1.0/dynrange)/activity['L'].max()) )
#    v2 = uint16( mri['t1_128']*(2**16/mri['t1_128'].max()) )

    #load Nifty data
    from nifti import NiftiImage
    v1 = NiftiImage('./activity_128.nii').data
    v2 = 0*v1

    #create volume renderer and initialize it
    V = VolumeRender((N,N,N),(512,512))
    V.show()
    V.set_volume1(v1)
    V.set_volume2(v2)

    #set visualization parameters
    V.set_density1(0.05)
    V.set_brightness1(7.1)
    V.set_transferOffset1(-0.04)
    V.set_transferScale1(1.1)
    V.set_density2(0.05)
    V.set_brightness2(0.3)
    V.set_transferOffset2(0.05)
    V.set_transferScale2(1)


    #visualize a dynamic scene and save screenshots
    N_frames = N*image_steps
    d_rotation_x = 0.0
    d_rotation_y = 0.5 
    d_zoom = 0.004

    im_frame_prev=0
    frame=0
    t_frame=0
    while 1:
        frame+=1
        if frame<1460:
            t_frame+=1
            if frame<150+250:
                V.rotate(d_rotation_x,d_rotation_y)
            if (frame>150+150 and frame<150+500):
                V.zoom(d_zoom)
            if frame>150+250:
                V.rotate(d_rotation_x/4,d_rotation_y/4)
            im_frame = t_frame/image_steps
            if not im_frame==im_frame_prev:
                im_frame_prev=im_frame
                if (im_frame<(N-thickness)):
                    v1b = double(v1)
                    v1b[:,:,im_frame:im_frame+thickness] = dynrange*double(v1[:,:,im_frame:im_frame+thickness])
                    #v1b = v1b*(2**16/v1b.max())
                    V.set_volume1(uint16(v1b))
        if frame>1500:
            t_frame-=1
            im_frame = t_frame/image_steps
            if not im_frame==im_frame_prev:
                im_frame_prev=im_frame
                if (im_frame<(N-thickness)):
                    v1b = double(v1)
                    v1b[:,:,im_frame:im_frame+thickness] = dynrange*double(v1[:,:,im_frame:im_frame+thickness])
                    #v1b = v1b*(2**16/v1b.max())
                    V.set_volume1(uint16(v1b))
        if (frame>2550+100):
            break
        elif (frame>2550+75):
            V.zoom(-d_zoom*3.5)
        elif frame>2550:
            V.rotate(d_rotation_x,d_rotation_y)
            V.zoom(-d_zoom*3.5)
        V.dump_screenshot("./screenshots/%d.png"%frame)
        sleep(0.003)

    while 1:
        pass


if __name__ == "__main__":
    run()





