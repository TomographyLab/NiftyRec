
from xmlrpclib import ServerProxy, Binary

DEF_ADDRESS = "http://localhost:10000"


class VolumeRender:
    def __init__(self,address=None):
        if address:
            self.server_address = address
        else:
            self.server_address = DEF_ADDRESS
        self.proxy = ServerProxy(self.server_address)

    def set_density1(self,density):
        return self.proxy.set_density1(density)

    def set_density2(self,density):
        return self.proxy.set_density2(density)

    def set_brightness1(self,brightness):
        return self.proxy.set_brightness1(brightness)

    def set_brightness2(self,brightness):
        return self.proxy.set_brightness2(brightness)

    def set_transferOffset1(self,transferOffset):
        return self.proxy.set_transferOffset1(transferOffset)

    def set_transferOffset2(self,transferOffset):
        return self.proxy.set_transferOffset2(transferOffset)

    def rotate(self,x,y):
        return self.proxy.rotate(x,y)

    def zoom(self,z):
        return self.proxy.zoom(z)

    def translate(self,x,y):
        return self.proxy.translate(x,y)

    def dump_screenshot(self,filename):
        self.proxy.dump_screenshot(filename)

    def set_volume1(self,volume):
        self.proxy.set_volume1(Binary(volume.tostring()))

    def set_volume2(self,volume):
        self.proxy.set_volume2(Binary(volume.tostring()))

    def set_colormap1(self,colormap):
        self.proxy.set_colormap1(Binary(colormap.tostring()))

    def set_colormap2(self,colormap):
        self.proxy.set_colormap2(Binary(colormap.tostring()))




