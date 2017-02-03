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
import AstroCell
plt.ioff()





class Image():
    """ Class that stores an image array, and provides various methods to process that image in some way;the output of these methods will be a new
    attribute of the Image object in question, each of which will *itslf* be a new Image object, with the same processing functionality """



    def __init__(self, in_image):
        """ Initialise the Image object """

        # Store the input image
        self.map = in_image



    def ScaleNorm(self):
        """ A function to normalise image scale so that the full dynamic range is taken advantage of """

        # Make copy of image to work with, and ensure array is float type
        out_image = self.map.copy()
        out_image = out_image.astype(float)

        # Subtract minimum pixel value from image, so that it becomes zero
        out_image -= np.nanmin(out_image)

        # Rescale map by a so that maximum pixel value is white (ie, value of 225)
        out_image *= 225.0 / np.nanmax(out_image)

        # Overwrite updated map with cleaned version
        self.map = out_image

        # Return updated Image
        return self





    def CleanEdges(self):
        """ Function  that removes rows/columns of single-value pixels from edge of an image. The 'unclean' image is saved to an attribute """

        # Make copy of image to work with
        in_image = self.map.copy()

        # Recrod Image object of uncleaned map as new attribute, for future reference
        self.unclean = Image(self.map)

        # Ensure that input image is int type
        in_image = np.round(in_image).astype(int)

        # Identify corner coordinates
        corner_coords = [[0, 0],
                        [0, in_image.shape[1]-1],
                        [in_image.shape[0]-1, 0],
                        [in_image.shape[0]-1, in_image.shape[1]-1]]

        # Create output image
        out_image = in_image.copy().astype(float)

        # Loop over corners
        for corner_coord in corner_coords:

            # Identify value of corner pixel, and create binary array for this value
            corner_pix = in_image[corner_coord[0], corner_coord[1]]
            bin_image = np.zeros(in_image.shape)
            bin_image[ np.where( in_image == corner_pix ) ] = 1

            # Identify any contiguous feature associated with corner pixel
            cont_structure = np.array([[0,1,0], [1,1,1], [0,1,0]])
            cont_image = scipy.ndimage.measurements.label(bin_image, structure=cont_structure)[0]
            cont_target = cont_image[corner_coord[0], corner_coord[1]]
            bin_image[ np.where( cont_image != cont_target ) ] = 0

            # Set pixels associated with any contiguous features to be NaN
            out_image[ np.where( cont_image == 1 ) ] = np.NaN

        # Overwrite existing map with cleaned version
        self.map = out_image

        # Return updated Image
        return self




















#class rgb():
#    """ Class for rgb image data and associated attributes """
#
#    def __init__(self, image_path):
#        """ Initialise rgb object """
#
#        # Read in the image, and convert into an array
#        bitmap_image = PIL.Image.open(image_path)
#        self.rgb = np.array(bitmap_image)
#
#        # Create a channel class for each band in the input array
#        bands = ['r','g','b']
#        for b in bands:
#            band = bands[b]
#            self.vars()[band] = channel()
#"""
#
#"""
#class channel():
#    """ Class for an image in a given band (ie, channel) """
#
#    def __init__(self, band, in_image):
#        """ Initialise band object """