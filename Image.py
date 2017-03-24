# Import smorgasbord
import pdb
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
from ChrisFuncs import SigmaClip
import AstroCell.Process
plt.ioff()





class RGB():
    """ Class that stores an RGB triplet of image arrays; the methods of an RGB object concern handling all three channels together"""



    def __init__(self, in_path):
        """ Initialise the RGB object """

        # Read in the image file, and convert to array
        bitmap_image = PIL.Image.open(in_path)
        rgb_image = np.array(bitmap_image)

        # Store the full rgb cube
        self.cube = rgb_image

        # Initialise AstroCell.Image object for each channel
        self.r = Image(rgb_image[:,:,0])
        self.g = Image(rgb_image[:,:,1])
        self.b = Image(rgb_image[:,:,2])

        # Create tuple containing each channel's Image object for iterating over
        self.iter = (self.r,self.g,self.b)

        # Make tuple of strings giving name of each channel, to use for iteration
        self.channels = ('r','g','b')



    def BlackOnWhite(self):
        """ A method that determines whether an image is on a black background; if not, image is inverted so that it is """

        # Iterate over rgb, performing rough canny extration
        cube_canny = np.zeros(self.cube.shape)
        for band in self.Iter():
            band.CannyCells()


#        # Get pixels values of image, with black (0) and white (255) pixels removed
#        in_values = self.cube.copy().flatten()
#        in_values = in_values[np.where(in_values!=0)]
#        in_values = in_values[np.where(in_values!=255)]
#        in_values = in_values[np.where(np.isnan(in_values)==False)]
#
#        r_canny = skimage.feature.canny(self.r.map, sigma=1.0)
#        r_canny_label = skimage.measure.label(np.invert(r_canny), connectivity=1)
#        astropy.io.fits.writeto('/home/chris/r_canny_label.fits', r_canny_label.astype(float), clobber=True)
#
#
#
#
#        # Get distribution of pixel values, and which values are over/under median
#        in_clip = SigmaClip(in_values, median=True, tolerance=0.0001, sigma_thresh=0.5)
#
#        # Construct histogram to find peak of pixel value distribution (bin size 10x Freedman-Diaconis bin size)
#        hist_bin_width = 10.0 * 2.0 * scipy.stats.iqr(in_values) * in_values.size**(-1./3.)
#        hist_bins = int( np.round( (np.nanmax(in_values)-np.nanmin(in_values)) / hist_bin_width ) )
#        pdb.set_trace()
#        hist = np.histogram(in_values, bins=hist_bins)
#
#
#        in_over = in_values[np.where(in_values>in_clip[1])]
#        in_under = in_values[np.where(in_values<in_clip[1])]
#        pdb.set_trace()





class Image():
    """ Class that stores an image array, and provides various methods to process that image in some way; the output of these methods can be a new
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




















