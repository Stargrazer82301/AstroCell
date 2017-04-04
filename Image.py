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
import scipy.ndimage.measurements
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



    def CannyCells(self, sigma=False):
        """ Wrapper around Process.CannyCells function """

        # Call Process.CannyCells function, with or without sigma kwarg, as approproate
        if sigma==False:
            canny_features = AstroCell.Process.CannyCells(self.map)
        elif isinstance(sigma, (int,float)):
            canny_features = AstroCell.Process.CannyCells(self.map, sigma=sigma)
        else:
            raise Exception('Sigma value to be passed to Canny filter is not neither a float nor an int.')

        # Record Canny feature map to Image object
        self.canny_features = canny_features




    def CannyCellStack(self):
        """ Method that stacks upon positions of identified features, to create a matched filter """

        # Identify segment indices present in map (excluding bakground feature)
        loop_features = list(set(self.canny_features.flatten()))[1:]

        # Calculate approximate diameters of features (assuming circular features)
        canny_areas = np.unique(self.canny_features, return_counts=True)[1].astype(float)
        canny_diams = 2.0 * np.sqrt(canny_areas/np.pi)

        # Hence decide size for cutout
        cutout_diam = 3.0*np.percentile(canny_diams,90.0)
        cutout_rad = int(np.round(0.5*cutout_diam))
        cutout_diam = (2*cutout_rad)+1

        # Create stack of NaNs to hold cutouts
        stack = np.NaN * np.zeros([int(cutout_diam), int(cutout_diam), len(loop_features)])

        # Create a padded version of the features map (and binary version of it), to enable creating cutouts
        pad_map = np.pad(self.map.astype(float), pad_width=cutout_rad+1, mode='constant', constant_values=np.NaN)
        pad_canny_features = np.pad(self.canny_features.astype(float), pad_width=cutout_rad+1, mode='constant', constant_values=np.NaN)
        pad_canny_bin = pad_canny_features.copy()
        pad_canny_bin[np.where(pad_canny_bin>0)] = 1


        # Loop over features
        for i in range(0,len(loop_features)):
            feature = loop_features[i]

            # Identify centre coords of feature
            feature_centre = scipy.ndimage.measurements.center_of_mass(pad_canny_bin, labels=pad_canny_features, index=feature)
            i_centre, j_centre = int(np.round(feature_centre[0])), int(np.round(feature_centre[1]))

            # Produce cutout centred on feature (rotate it an arbitary number of times?), and insert into stack
            cutout = pad_map[i_centre-cutout_rad:i_centre+cutout_rad+1, j_centre-cutout_rad:j_centre+cutout_rad+1]
            #cutout = np.rot90( cutout, k=np.random.randint(0, high=4) )
            stack[:,:,i] =cutout

        # Make stack coadd, and record to object
        stack_coadd = np.nanmean(stack, axis=2)
        #astropy.io.fits.writeto('/home/chris/stack_coadd.fits', stack_coadd, clobber=True)
        self.canny_stack = stack_coadd







