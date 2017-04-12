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
        kernel = astropy.convolution.kernels.Gaussian2DKernel(5.0)
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
        pad_map = np.pad(self.detmap.astype(float), pad_width=cutout_rad+1, mode='constant', constant_values=np.NaN)
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
            stack[:,:,i] = cutout

        # Make stack coadd
        stack_coadd = np.nanmean(stack, axis=2)

        # Normalise stack, and make edges equal to zero (so it works better as a kernel)
        stack_edges = np.array([stack_coadd[0,:], stack_coadd[-1,:], stack_coadd[:,0], stack_coadd[:,-1]]).flatten()
        stack_coadd -= np.nanmedian(stack_edges)
        stack_coadd /= np.nansum(stack_coadd)

        # Record stack
        self.canny_stack = stack_coadd
        #astropy.io.fits.writeto('/home/chris/stack_coadd.fits', stack_coadd, clobber=True)



    def ThreshSegment(self, bg_mask=None):
        """ A method that uses basic threshold segmenting to identify cells """

        # Use areas of this channel's Canny features to decide minimum pixel area limit for segments
        canny_areas = np.unique(self.canny_features, return_counts=True)[1].astype(float)
        canny_areas_clip = SigmaClip(canny_areas, median=True, sigma_thresh=2.0)
        area_thresh = int( np.round( canny_areas_clip[1] - ( 2.0 * canny_areas_clip[0] ) ) )
        #canny_diam = 2.0 * np.sqrt(area_thresh/np.pi)

        # If no features smaller than peak (ie, the modal size is also the smallest size), set this value to be the threshold
        if np.isnan(area_thresh):
            area_thresh = canny_areas_clip[1]

        # Perform sigma clipping, to charactarise background
        in_map = self.detmap.copy().astype(float)
        bg_map = in_map.copy()
        if bg_mask!=None:
            bg_map[ np.where(bg_mask>0) ] = np.NaN
        bg_clip = SigmaClip(bg_map, median=False, sigma_thresh=3.0)

        # Background subtract map, and determine segmentation threshold
        in_map -= bg_clip[1]
        #seg_thresh = skimage.filters.threshold_otsu(in_map, nbins=1024)
        seg_thresh = 1.5 * bg_clip[0]

        # Use photutils to segment map
        seg_map = photutils.detect_sources(in_map, threshold=seg_thresh, npixels=area_thresh, connectivity=8).array
        seg_map = AstroCell.Process.LabelShuffle(seg_map)

        """# Crude test of deblending
        seg_areas = np.unique(seg_map, return_counts=True)[1].astype(float)
        conv_map = astropy.convolution.convolve_fft(self.r.detmap.copy(), rgb.r.canny_stack, interpolate_nan=True, normalize_kernel=True, boundary='reflect', allow_huge=True)
        tophat_kernel = astropy.convolution.Tophat2DKernel( (np.median(canny_areas)/np.pi)**0.5 )
        deblend_map = photutils.deblend_sources(in_map, seg_map, npixels=area_thresh, filter_kernel=tophat_kernel, labels=None, nlevels=1024, contrast=0.000000001, mode='exponential', connectivity=8, relabel=True)
        grey_erode_structure = scipy.ndimage.generate_binary_structure(2,1)
        grey_erode_map = scipy.ndimage.grey_erosion(in_map, structure=grey_erode_structure, mode='reflect')
        astropy.io.fits.writeto('/home/chris/red_grey_erode.fits', grey_erode_map, clobber=True)"""

        # Record attributes
        self.thresh_segmap = seg_map
        self.thresh_level = seg_thresh
        self.thresh_area = area_thresh






