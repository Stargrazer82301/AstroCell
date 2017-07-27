# Import smorgasbord
import pdb
import sys
import os
import shutil
import warnings
import multiprocessing as mp
import numpy as np
import scipy.stats
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.io.fits
import skimage.segmentation
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



class Parallel():
    """ Class that stores whether, and what kind, of parallel operation is desired """

    def __init__(self, parallel):
        """ Initialise the Parallel object """

        # Parse and record parallel request information
        if parallel == False:
            self.parallel = False
            self.threads = 1
            self.subthreads = 1
        elif isinstance(parallel,float) or isinstance(parallel,float):
            self.parallel = True
            self.threads = 1
            self.subthreads = int(parallel)
        elif parallel == True:
            self.parallel = True
            self.threads = 1
            self.subthreads = mp.cpu_count()-1




    # Plot label borders onto image, for each label
    #astropy.io.fits.writeto('/home/chris/bob.fits', image_label_dilate.astype(float), clobber=True)



class Bunch:
    """ Convenience class for gathering related data in an object """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


