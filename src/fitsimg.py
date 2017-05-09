import glob
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch, SqrtStretch, AsinhStretch, HistEqStretch,ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt


class FitsImage(object):
    """
    A helper class when reading fits files. Depends on pyFits.

    """
    DIRLOC = ''
    
    def __init__(self,filename):
        self.filename = filename
        self.hdulist = fits.open(self.filename)
        #print self.hdulist.info()
        self.header = self.hdulist[0].header
        data = self.hdulist[0].data
        self.data = data.astype(float)
        #print self.data.shape
    
    def resize(self, shape):
        """
        A function to resize the image data to a shape, preserving the range.
        
        INPUT:
            shape - tuple
        
        EXAMPLE:
            f = FitsImage("fits_10_blue_red_ir_1500/0001_266.6147_27.72068_dss2blue_N333.fits")
            f.resize((100,100))
            f.plot()
            f.data.shape
        """
        self.data = resize(self.data,shape,preserve_range=True)
        
    def cropcenter(self,w,h,points=False):
        """
        Returns an image array around the center of an image array.
        """
        shape = self.data.shape
        center = (shape[0] / 2, shape[1] / 2)
        hl = center[1]-h/2
        hr = center[1]+h/2
        wl = center[0]-w/2
        wr = center[0]+w/2
        self.data = self.data[hl:hr,wl:wr]
        #if points == False:
        #    return self.data_crop
        #else:
        #    return self.data_crop, (wl,wl,wr,wr), (hl,hr,hl,hr)
        
    def plot(self,stretch="hist",cmap="gray",origin="lower"):
        self.fig, self.ax = plt.subplots()
        if stretch=="hist":
            norm = ImageNormalize(stretch=HistEqStretch(self.data))
            self.im = self.ax.imshow(self.data,cmap=cmap,origin=origin,norm=norm)
        else:
            self.im = self.ax.imshow(self.data,cmap=cmap,origin=origin)
        self.ax.set_xlim(0,self.data.shape[0])
        self.ax.set_ylim(0,self.data.shape[1])
        self.fig.colorbar(self.im)
        
    def hist(self,num=1000,figsize=(20,8)):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=1,figsize=figsize)
        self.axes.hist(self.data.flat,num)
        
    def plot_hist(self,cmap='gray',figsize=(20,8),vmin=None,vmax=None,scaled=False,num=1000,titlepre=""):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=2,figsize=figsize)
        im = self.axes[0].imshow(self.data,cmap=cmap,vmin=vmin,vmax=vmax)
        #self.axes[0].set_title(titlepre + "DATE:" + self.header["DATE-OBS"] + " EXP: " + str(self.header["EXPOSURE"]) + " SHAPE: " +str(self.data.shape))
        self.fig.colorbar(im)
        self.axes[1].hist(self.data.flat,num,log=True)
        
    def close(self):
        self.hdulist.close()
