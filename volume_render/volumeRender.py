
from ctypes import CDLL, c_float, c_uint, pointer, POINTER, c_uint8, c_int, c_ubyte
from thread import start_new_thread
from numpy import float, float32, ascontiguousarray, uint16, uint8, int32, uint, fromstring
from Image import frombuffer
from time import sleep
from SimpleXMLRPCServer import SimpleXMLRPCServer
import platform 

DEF_PORT  = 10000

if platform.system() == 'Linux':
    LIB_NAME = 'libvolumeRender.so'
elif platform.system() == 'Darwin':
    LIB_NAME = 'libvolumeRender.dylib'
elif platform.system() == 'Windows':
    LIB_NAME = 'volumeRender.dll' 

class VolumeRender:
    def __init__(self,volume_size,camera_size,dtype=uint16,port=DEF_PORT):
        #try:
        self.lib = CDLL(LIB_NAME)
        #except OSError:
        #    self.lib=None
        self.volume_size = volume_size
        self.camera_size = camera_size
        self.gui_thread=None
        #Check if requested dtype is a compatible type
        if dtype not in [uint8, uint16, float32]:
            print "Can't handle the specified dtype: ",str(dtype)
            self.lib=None
            self.dtype=None
        self.dtype=dtype
        self.port=port
        sleep(0.3)

    def start_server(self):
        self.server = SimpleXMLRPCServer(('',self.port))
        self.server.register_function(self.set_density1)
        self.server.register_function(self.set_density2)
        self.server.register_function(self.set_brightness1)
        self.server.register_function(self.set_brightness2)
        self.server.register_function(self.set_volume1_xmlrpc, 'set_volume1')
        self.server.register_function(self.set_volume2_xmlrpc, 'set_volume2')
        self.server.register_function(self.set_colormap1_xmlrpc, 'set_colormap1')
        self.server.register_function(self.set_colormap2_xmlrpc, 'set_colormap2')
        self.server.register_function(self.set_transferOffset1)
        self.server.register_function(self.set_transferOffset2)
        self.server.register_function(self.rotate)
        self.server.register_function(self.translate)
        self.server.register_function(self.zoom)
        self.server.register_function(self.dump_screenshot)
        start_new_thread(self.server.serve_forever,())
        
    def show(self,as_server=False):
        if self.lib==None:
            print "Couldn't load library"
            return True
        params = ( int(self.volume_size[0]),int(self.volume_size[1]),int(self.volume_size[2]),int(self.camera_size[0]),int(self.camera_size[1]) )
        if as_server:
            self.start_server()
            self.lib.run_gui( int(self.volume_size[0]),int(self.volume_size[1]),int(self.volume_size[2]),int(self.camera_size[0]),int(self.camera_size[1]) )
        else:
            self.gui_thread = start_new_thread(self.lib.run_gui,params)
        sleep(0.3)
        return False

    def get_density1(self):
        return self.lib.get_density1()

    def set_density1(self,density):
        return self.lib.set_density1(c_float(density))

    def get_brightness1(self):
        return self.lib.get_brightness1()

    def set_brightness1(self,brightness):
        return self.lib.set_brightness1(c_float(brightness))

    def get_transferOffset1(self):
        return self.lib.get_transferOffset1()

    def set_transferOffset1(self,transferOffset):
        return self.lib.set_transferOffset1(c_float(transferOffset))

    def get_transferScale1(self):
        return self.lib.get_transferScale1()

    def set_transferScale1(self,transferScale):
        return self.lib.set_transferScale1(c_float(transferScale))        


    def get_density2(self):
        return self.lib.get_density2()

    def set_density2(self,density):
        return self.lib.set_density2(c_float(density))

    def get_brightness2(self):
        return self.lib.get_brightness2()

    def set_brightness2(self,brightness):
        return self.lib.set_brightness2(c_float(brightness))

    def get_transferOffset2(self):
        return self.lib.get_transferOffset2()

    def set_transferOffset2(self,transferOffset):
        return self.lib.set_transferOffset2(c_float(transferOffset))

    def get_transferScale2(self):
        return self.lib.get_transferScale2()

    def set_transferScale2(self,transferScale):
        return self.lib.set_transferScale2(c_float(transferScale))    


    def set_volume1(self,volume):        
        #Convert volume to self.dtype
        #if not(volume.dtype==self.dtype):
        volume = ascontiguousarray(volume,dtype=self.dtype)
        print volume.ctypes.data
        self.lib.set_volume1.argtypes = [POINTER(c_ubyte)]
        return self.lib.set_volume1(volume.ctypes.data_as(POINTER(c_ubyte)))

    def set_volume2(self,volume):        
        #Convert volume to self.dtype
        #if not(volume.dtype==self.dtype):
        volume = ascontiguousarray(volume,dtype=self.dtype)
        self.lib.set_volume2.argtypes = [POINTER(c_ubyte)]
        return self.lib.set_volume2(volume.ctypes.data_as(POINTER(c_ubyte)))

    def set_volume1_xmlrpc(self,volume):
        return self.set_volume1(fromstring(volume.data,dtype=self.dtype))

    def set_volume2_xmlrpc(self,volume):
        return self.set_volume2(fromstring(volume.data,dtype=self.dtype))

    def set_colormap1(self,colormap):
        #Convert volume to float32
        if not(colormap.dtype==float32):
            colormap = ascontiguousarray(colormap,dtype=float32)
        if not colormap.shape[1]==4:
            print "Colormap must be of size Nx4"
        colormap_size = uint(colormap.shape[0])
        self.lib.set_colormap1.argtypes = [POINTER(c_float), c_uint]
        self.lib.set_colormap1.restype = c_int
        colormap_p = colormap.ctypes.data_as(POINTER(c_float))
        self.lib.set_colormap1(colormap_p,colormap_size)
        return True

    def set_colormap2(self,colormap):
        #Convert volume to float32
        if not(colormap.dtype==float32):
            colormap = ascontiguousarray(colormap,dtype=float32)
        if not colormap.shape[1]==4:
            print "Colormap must be of size Nx4"
        colormap_size = uint(colormap.shape[0])
        self.lib.set_colormap2.argtypes = [POINTER(c_float), c_uint]
        self.lib.set_colormap2.restype = c_int
        colormap_p = colormap.ctypes.data_as(POINTER(c_float))
        self.lib.set_colormap2(colormap_p,colormap_size)
        return True

    def set_colormap1_xmlrpc(self,colormap):
        colormap = fromstring(colormap.data)
        colormap = colormap.reshape(colormap.size/4,4)
        return self.set_colormap1(colormap)

    def set_colormap2_xmlrpc(self,colormap):
        colormap = fromstring(colormap.data)
        colormap = colormap.reshape(colormap.size/4,4)
        return self.set_colormap2(colormap)

    def set_raystep(self,step):
        return self.lib.setRayStep(c_float(step))

#    def get_screenshot(self):
#        x = c_uint(); y = c_uint();
#        p_x = pointer(x); p_y = pointer(y);
#        self.lib.get_screenshot.argtypes = [POINTER(c_uint),POINTER(c_uint)]
#        self.lib.get_screenshot.restype = POINTER(c_uint8)
#        p_image = self.lib.get_screenshot(p_x, p_y);
#        self.lib.get_screenshot.restype = POINTER(c_uint8*3*int32(x)*int32(y))
#        p_image = self.lib.get_screenshot(p_x, p_y);
        
    def dump_screenshot(self,filename):
        self.lib.dump_screenshot()
        fid = open('./screenshot.raw')
        buff = fid.read()
        fid.close()
        im = frombuffer('RGB',(self.camera_size[0],self.camera_size[1]),buff)
        im.save(filename)
        return True

    def translate(self,x,y):
        self.lib.translate.argtypes = [c_float,c_float]
        self.lib.translate(float(x),float(y))
        return True
		
    def rotate(self,x,y):
        self.lib.rotate.argtypes = [c_float,c_float]
        self.lib.rotate(float(x),float(y))
        return True

    def zoom(self,z):
        self.lib.zoom.argtypes = [c_float]
        self.lib.zoom(float(z))
        return True
		
    def __del__(self):
        self.stop()

    def stop(self):
        self.lib.stop()
        if self.gui_thread:
            self.gui_thread.join()



if __name__ == "__main__":
    import sys
    volume_size = (128,128,128)
    camera_size = (512,512)
    if len(sys.argv)>=4:
        volume_size = (int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
    if len(sys.argv)>=6:
        camera_size = (int(sys.argv[4]),int(sys.argv[5]))        
    renderer = VolumeRender(volume_size, camera_size, uint16, DEF_PORT)
    renderer.show(as_server=True)



