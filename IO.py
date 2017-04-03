# Import smorgasbord
import pdb
import sys
import os
import shutil
import numpy as np
import scipy.stats
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.io.fits
import PIL.Image
import AstroCell.Image
plt.ioff()





def LoadRGB(in_path):
    """ Function that reads in an image file and returns three AstrCell.Image objects (r, g, b) """

    # Read in image, and conver to array
    bitmap_image = PIL.Image.open(in_path)
    rgb_image = np.array(bitmap_image)

    # Create an Image object from each channel in turn
    r = AstroCell.Image.Image(rgb_image[:,:,0])
    g = AstroCell.Image.Image(rgb_image[:,:,1])
    b = AstroCell.Image.Image(rgb_image[:,:,2])

    # Return Image objects
    return r, g, b



class Bunch:
    """ Convenience class for gathering related data """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)



class TempDir():
    """ Convenience class for working with temporary directory"""

    def __init__(self, out_dir):
        """ Initialise the temporary directory object """
        self.dir = os.path.join(out_dir,'Temp')
        if not os.path.exists(self.dir):
            os.mkdir(self.dir)

    def qw(self, data, name):
        """ Method for quickly writing a numpy array to a FITS file in the temporary directory """
        if name[-4:] != '.fits':
            if name[-1:] == '.':
                name = name[:-1] + '.fits'
            else:
                name += '.fits'
        if data.dtype == 'bool':
            data = data.astype(int)
        astropy.io.fits.writeto(os.path.join(self.dir,name), data, clobber=True)


