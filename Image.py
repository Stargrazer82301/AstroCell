# Import smorgasbord
import pdb
import time
import copy
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
from ChrisFuncs.Photom import EllipseMask
import AstroCell.Process
plt.ioff()





class Image():
    """ Class that stores an image array, and provides various methods to process that image in some way; the output of these methods can be a new
    attribute of the Image object in question, each of which will *itslf* be a new Image object, with the same processing functionality """



    def __init__(self, in_image):
        """ Initialise the Image object """

        # Store the input image
        self.map = in_image.astype(float)



    def CleanEdges(self):
        """ Method that removes rows/columns of single-value pixels from edge of an image. The 'unclean' image is saved to an attribute """

        # Make copy of image to work with
        in_image = self.map.copy()

        # Recrod Image object of uncleaned map as new attribute, for future reference
        self.unclean = Image(self.map)

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

        # Loop over rows, looking for rows which all contain the same value
        for i in range(0,in_image.shape[0]):
            if np.nanstd(in_image[i,:]) == 0:
                out_image[i,:] = np.array([np.nan]*in_image.shape[1])

        # Loop over columns, looking for rows which all contain the same value
        for i in range(0,in_image.shape[1]):
            if np.nanstd(in_image[:,i]) == 0:
                out_image[:,i] = np.array([np.nan]*in_image.shape[0])

        # Use a convolved version of the input map, with interpolation over NaN pixels, to impute replacement values for NaNs
        kernel = astropy.convolution.kernels.Gaussian2DKernel(1.5)
        conv_map = astropy.convolution.convolve_fft(out_image, kernel, interpolate_nan=True, normalize_kernel=True, boundary='reflect')
        out_image[np.where(np.isnan(out_image))] = conv_map[np.where(np.isnan(out_image))]

        # Overwrite existing map with cleaned version
        self.map = out_image





class RGB():
    """ Class that stores an RGB triplet of image arrays; the methods of an RGB object concern handling all three channels together"""



    def __init__(self, in_path):
        """ Initialise the RGB object """

        # Store path for future use
        self.path = in_path

        # Read in the image file, and convert to array
        bitmap_image = PIL.Image.open(in_path)
        rgb_image = np.array(bitmap_image)

        # Store raw, unaltered versionf of each channel as AstroCell.Image objects
        class Raw(object):
            pass
        self.raw = Raw()
        self.raw.r = Image(rgb_image[:,:,0].copy())
        self.raw.g = Image(rgb_image[:,:,1].copy())
        self.raw.b = Image(rgb_image[:,:,2].copy())

        # Store the full rgb cube, as a float
        self.cube = rgb_image.astype(float)

        # Initialise primary AstroCell.Image object for each channel
        self.r = Image(self.cube[:,:,0].copy())
        self.g = Image(self.cube[:,:,1].copy())
        self.b = Image(self.cube[:,:,2].copy())

        # Create tuple containing each channel's Image object for iterating over
        self.iter = (self.r,self.g,self.b)

        # Make tuple of strings giving name of each channel, to use for iteration
        self.channels = ('r','g','b')

        # Set data type for each channel to be float
        for channel in self.iter:
            channel.map = channel.map.astype(float)

        # Record each channel's name as an attribute
        for i in range(0, len(self.iter)):
            self.iter[i].name =  self.channels[i]



    def MakeCoadd(self):
        """ Method that adds an new Image object, that's a coadd of the three channels """

        # Initialise AstroCell.Image object for a coadd of all three channels
        coadd = np.sum(self.cube,axis=2).astype(float) / 3.0
        self.coadd = Image(coadd)

        # Create tuple containing each channel's Image object (including the coadd) for iterating over
        self.iter_coadd = (self.r,self.g,self.b,self.coadd)

        # Make tuple of strings giving name of each channel (including the coadd), to use for iteration
        self.channels = ('r','g','b','coadd')

        # Record name of coadd channel
        self.coadd.name = 'coadd'



    def BlackOnWhite(self):
        """ A method that determines whether an image is on a black background; if not, image is inverted so that it is """

        # Do a rough Canny-based cell extraction on each channel, and coadd the Canny feature maps
        canny_cube = np.zeros(self.cube.shape)
        for i in range(0, len(self.iter)):
            canny_cube[:,:,i] = AstroCell.Process.CannyCells(self.iter[i].map)
        self.canny_features = np.sum(canny_cube, axis=2)
        canny_where = np.where(self.canny_features>0)
        #astropy.io.fits.writeto('/home/chris/canny_features.fits', self.canny_features, clobber=True)

        # Flag pixel values for pixels within Canny features
        cube = self.cube.copy().astype(int)
        canny_where = np.where(self.canny_features>0)
        for i in range(0,cube.shape[2]):
            cube[:,:,i][canny_where] = -99

        # Get pixel values of image, with Canny-flagged (-99), black (0), and white (255) pixels removed
        values = cube.copy().flatten().astype(float)
        values = values[np.where(values!=-99)]
        values = values[np.where(values>0)]
        values = values[np.where(values<255)]
        values = values[np.where(np.isnan(values)==False)]

        # Construct histogram to find peak of pixel value distribution (using Freedman-Diaconis rule to decide bin size)
        hist_bin_width = 2.0 * scipy.stats.iqr(values) * values.size**(-1./3.)
        hist_bins = int( np.round( (np.nanmax(values)-np.nanmin(values)) / hist_bin_width ) )
        hist = np.histogram(values, bins=hist_bins)
        hist_peak = 0.5 * ( hist[1][np.argmax(hist[0])] + hist[1][np.argmax(hist[0])+1] )

        # Find number of values above and below peak of pixel value distribution
        values_above = values[ np.where(values>hist_peak) ]
        values_below = values[ np.where(values<hist_peak) ]

        # If distribution has strong negative skew, assume the image has black background
        skew = len(values_below)/len(values_above)
        if skew<0.5:
            self.inverted = False

        # Else assume image has white background, and invert Image
        else:
            self.inverted = True
            for channel in self.iter:
                channel.map = -1.0 * ( channel.map - 255.0 )
            self.coadd.map = -1.0 * ( self.coadd.map - 255.0 )




    def DetFilter(self):
        """ Method that removes smooth all large-scale background structure from map, to create an optimal detection map """

        # Use canny features map to get distribution of filter sizes
        canny_diams = []
        canny_features_flat = self.canny_features.flatten()
        for i in range(1, int(np.max(canny_features_flat))):
            canny_diam = 2.0 * np.sqrt( float(np.where(canny_features_flat==i)[0].shape[0]) / np.pi)
            canny_diams.append(canny_diam)
        canny_diams = np.array(canny_diams)[np.where(np.array(canny_diams)>0)]

        # Decide size of of minimum filter to apply, and shape of filter footprint
        kernel_size = 4.0 * np.percentile(canny_diams, 90.0)
        kernel = astropy.convolution.kernels.Gaussian2DKernel(kernel_size)



        canny_bg_map = self.coadd.map.copy().astype(float)
        canny_bg_map[ np.where(self.canny_features>0) ] = np.NaN
        conv_map = astropy.convolution.convolve_fft(canny_bg_map, kernel, interpolate_nan=True, normalize_kernel=True, boundary='reflect', allow_huge=True)

        #pdb.set_trace()
        conv_map[ np.where( np.isnan(self.r.map)==True ) ] = np.NaN
        conv_sub = self.r.map - conv_map
        #astropy.io.fits.writeto('/home/chris/conv.fits', conv_map, clobber=True)

        # Iterate over each channel, applying a minimum filter to each
        for channel in self.iter:
            pdb.set_trace()














