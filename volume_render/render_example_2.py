
from volumeRender import VolumeRender
from numpy import ones, zeros, uint16, double
from time import sleep
import scipy.io
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
    V.set_brightness1(5.4)
    V.set_transferOffset1(-0.02)
    V.set_transferScale1(1.27)
    V.set_density2(0.46)
    V.set_brightness2(0.5)
    V.set_transferOffset2(0.06)
    V.set_transferScale2(1.31)

    sleep(10)

    for im_frame in range(126):
        v1b = double(v1)
        v1b[:,im_frame:im_frame+thickness,:] = dynrange*double(v1[:,im_frame:im_frame+thickness,:])
        V.set_volume1(uint16(v1b))
        V.dump_screenshot("./screenshots/example2_%d.png"%im_frame)

    while 1:
        pass



if __name__ == "__main__":
    run()

