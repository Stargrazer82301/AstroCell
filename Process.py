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



def EdgeClean(in_image):
    """ Function that removes rows/columns of single-value pixels from edge of an image """

    # Ensure that input image is in int format
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

    # Return corrected image to user
    return out_image


