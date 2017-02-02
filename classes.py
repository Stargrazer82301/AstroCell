# Import smorgasbord
import pdb
import sys
import os
import shutil
import imghdr
import inspect
import numpy as np
import scipy.stats
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.convolution
import astropy.stats
import astropy.visualization
import astropy.visualization.mpl_normalize
import astropy.io.fits
import photutils
import skimage.feature
import PIL.Image
import GalCell
plt.ioff()



#class rgb():
#    """ Holder class for rgb image data and associated attributes """
#
#    def __init__(self, image_path):







class temp_dir():
    """ Convenience class for working with temporary directory"""

    def __init__(self, out_dir):
        """ Initialise the temporary directory object """
        self.dir = os.path.join(out_dir,'Temp')
        if not os.path.exists(self.dir):
            os.mkdir(self.dir)

    def qw(self, data, name):
        """ Method for quickly writing a numpy array to a FITS file in the temporary directory """
        if name[-4:]!='.fits':
            if name[-1:]=='.':
                name = name[:-1] + '.fits'
            else:
                name += '.fits'
        astropy.io.fits.writeto(os.path.join(self.dir,name), data, clobber=True)

