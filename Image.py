# Import smorgasbord
import pdb
import copy
import gc
import multiprocessing as mp
import numpy as np
import scipy.stats
import scipy.ndimage.measurements
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
import skimage.segmentation
import joblib
from ChrisFuncs import SigmaClip
from ChrisFuncs.Photom import EllipseMask
import AstroCell.Process
plt.ioff()





class Image():
    """ Class that stores an image array, and provides various methods to process that image in some way """



    def __init__(self, in_image):
        """ Initialise the Image object """

        # Store the input image
        self.map = in_image.astype(float)



    def Raw(self):
        """ Short method that makes an untouched copy of the Image, nested under itself """

        # Initialise raw Image
        self.raw = self.map.copy()



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

        """# Use a convolved version of the input map, with interpolation over NaN pixels, to impute replacement values for NaNs
        kernel = astropy.convolution.kernels.Gaussian2DKernel(5.0)
        conv_map = astropy.convolution.convolve_fft(out_image, kernel, interpolate_nan=True, normalize_kernel=True, boundary='reflect')
        out_image[np.where(np.isnan(out_image))] = conv_map[np.where(np.isnan(out_image))]"""

        # Use map median to impute replacement value for NaNs
        out_image[np.where(np.isnan(out_image))] = np.nanmedian(in_image)

        # Overwrite existing map with cleaned version
        self.map = out_image



    def CannyBlobs(self, sigma=False):
        """ A method that uses Canny edge detection to generate a initial guess of what regions of an image are occupied by cells """

        # Create coy of map to work with
        in_map = self.map.copy()

        # Evaluate noise in image, via iterative sigma-clipping
        clip = SigmaClip(self.map, median=True, sigma_thresh=3.0)

        # Use noise level to generate atificial noise to add to image
        canny_iter = 100
        noise = np.random.normal(0.0, 0.25*clip[0], size=(in_map.shape[0],in_map.shape[1],canny_iter))

        # Add each generation of random noise to image in turn, performing Canny edge filtering on each (in parallel with joblist)
        canny_stack = np.zeros(noise.shape)
        if self.parallel:
            canny_list = joblib.Parallel( n_jobs=mp.cpu_count()-1 )\
                                        ( joblib.delayed(skimage.feature.canny)\
                                        (in_map+noise[:,:,i], sigma=sigma) \
                                        for i in range(0,noise.shape[2]) )
            for i in range(0,noise.shape[2]):
                canny_stack[:,:,i] = canny_list[i]
            del(canny_list)
        else:
            for i in range(0,noise.shape[2]):
                canny_stack[:,:,i] = skimage.feature.canny(in_map+noise[:,:,i], sigma=sigma)

        # Co-add each individual Canny iteration; Canny edges found in lots of iterations are assumed to be real
        canny_iter_frac = 0.35
        canny = np.sum(canny_stack, axis=2)
        canny[ np.where( canny < (canny_iter*canny_iter_frac) ) ] = 0
        canny[ np.where( canny >= (canny_iter*canny_iter_frac) ) ] = 1
        #astropy.io.fits.writeto('/home/chris/canny.fits', canny.astype(float), clobber=True)

        # Run map through one iteration of binary closing, to close up gaps in the perimeters of canny edges
        canny_close = scipy.ndimage.binary_closing(canny, iterations=1, structure=scipy.ndimage.generate_binary_structure(2,1))

        # Label closed Canny image features
        canny_cell_labels = skimage.measure.label(np.invert(canny_close), connectivity=1)

        # Record number of pixels in each labelled feature
        canny_cell_areas = np.unique(canny_cell_labels, return_counts=True)[1]
        canny_cell_areas = canny_cell_areas[np.where(canny_cell_areas>0)]

        # Clip label areas, to identify threshold above which cell regions are large enough to likely be spurious
        labels_clip = SigmaClip(canny_cell_areas, median=True, sigma_thresh=5.0)
        labels_area_thresh = np.max(labels_clip[2])
        labels_exclude = np.arange(0,canny_cell_areas.size)[ np.where( (canny_cell_areas>labels_area_thresh) | (canny_cell_areas<=5) ) ]

        # Remove spurious labels (flattening to improve speed, then reshaping after processing)
        canny_cells = canny_cell_labels.copy().flatten()
        canny_cells[np.in1d(canny_cells,labels_exclude)] = 0
        canny_cells = np.reshape(canny_cells, canny_cell_labels.shape)

        # Fill cell regions, dilate them by 1 pixel (to account for the width of the Canny border), then relabel
        canny_fill = scipy.ndimage.binary_fill_holes(canny_cells, structure=scipy.ndimage.generate_binary_structure(2,1))
        canny_dilate = scipy.ndimage.binary_dilation(canny_fill, iterations=1, structure=scipy.ndimage.generate_binary_structure(2,2))
        canny_cells = skimage.measure.label(canny_dilate, connectivity=1)
        canny_cells = AstroCell.Process.LabelShuffle(canny_cells)
        #astropy.io.fits.writeto('/home/chris/canny_cells.fits', canny_cells, clobber=True)

        # Catch and fix when most of map is a 'feature'
        mode = scipy.stats.mode(canny_cells.flatten())[0][0]
        if mode != 0:
            mode_frac = np.where(canny_cells==mode)[0].shape[0] / canny_cells.size
            if mode_frac > 0.5:
                canny_cells = np.zeros(canny_cells.shape).astype(int)

        # Record final image
        gc.collect()
        self.canny_features = canny_cells



    def CannyCellStack(self):
        """ Method (redundant?) that stacks upon positions of identified features, to create a matched filter """

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



    def LogDogBlobs(self, canny_features=None, force_attribute=False):
        """ A method that uses Laplacian-of-Gaussian and Difference-of-Gaussian blob detection to identify which pixels have cells in """

        # If canny features map provided, use this; otherwise just use features map for this channel
        if not isinstance(canny_features,np.ndarray):
            canny_features = self.canny_features

        # Use canny features in this channel to work out range of cell sizes
        if np.nanstd(canny_features) == 0:
            canny_areas = np.array([])
        else:
            canny_areas = np.unique(canny_features, return_counts=True)[1].astype(float)
        canny_diams = 2.0 * np.sqrt(canny_areas/np.pi)
        diams_clip = SigmaClip(canny_diams, median=True, sigma_thresh=3.0)
        diams_thresh_min = 0.5 * np.min(diams_clip[2])
        diams_thresh_max = 0.5 * (diams_clip[1] + diams_clip[0])

        # Run LoG extraction
        log_blobs = skimage.feature.blob_log(self.map.copy(), min_sigma=diams_thresh_min, max_sigma=diams_thresh_max,
                                             num_sigma=25, overlap=0.95, threshold=0.1)

        # Run DoG extraction
        dog_blobs = skimage.feature.blob_dog(self.map.copy(), min_sigma=diams_thresh_min, max_sigma=diams_thresh_max,
                                             sigma_ratio=1.3, overlap=0.95, threshold=0.1)

        # Convert third column of blob outputs to radii
        log_blobs[:,2] = log_blobs[:,2] * np.sqrt(2.0)
        dog_blobs[:,2] = dog_blobs[:,2] * np.sqrt(2.0)

        # Create mask
        blob_mask = np.zeros(self.map.shape)

        # Loop over LoG blobs, adding them to mask
        for i in range(0, log_blobs.shape[0]):
            blob_mask += EllipseMask(blob_mask, log_blobs[i,2], 1.0, 0.0, log_blobs[i,0], log_blobs[i,1])

        # Loop over DoG blobs, adding them to mask
        for i in range(0, dog_blobs.shape[0]):
            blob_mask += EllipseMask(blob_mask, dog_blobs[i,2], 1.0, 0.0, dog_blobs[i,0], dog_blobs[i,1])

        # Use mask to create a features map
        blob_features = scipy.ndimage.measurements.label(blob_mask)[0]
        blob_features = AstroCell.Process.LabelShuffle(blob_features)

        # Return mask
        if force_attribute:
            self.logdog_mask = blob_mask
            self.logdog_features = blob_features
        else:
            return blob_mask, blob_features
        #astropy.io.fits.writeto('/home/chris/blob_mask.fits', blob_mask, clobber=True)



    def ThreshSegment(self, bg_mask=None):
        """ A method that uses basic threshold segmenting to identify cells """

        # Use areas of this channel's Canny features to decide minimum pixel area limit for segments
        if np.nanstd(self.canny_features) == 0:
            area_thresh = 0.0
        else:
            canny_areas = np.unique(self.canny_features, return_counts=True)[1].astype(float)
            canny_areas_clip = SigmaClip(canny_areas, median=True, sigma_thresh=2.0)
            area_thresh = int( np.round( canny_areas_clip[1] - ( 3.0 * np.nanstd(canny_areas[np.where(canny_areas<canny_areas_clip[1])]) ) ) )
            #canny_diam = 2.0 * np.sqrt(area_thresh/np.pi)

        # If no features smaller than peak (ie, the modal size is also the smallest size), set this value to be the threshold
        if np.isnan(area_thresh):
            area_thresh = canny_areas_clip[1]

        # Handle maps with no features
        area_thresh = np.max([1.0, area_thresh])

        # Perform sigma clipping, to charactarise background
        bg_map = self.detmap.copy().astype(float)
        if isinstance(bg_mask,np.ndarray):
            bg_map[ np.where(bg_mask>0) ] = np.NaN
        bg_clip = SigmaClip(bg_map, median=False, sigma_thresh=2.0)

        # Background subtract map, and determine segmentation threshold
        in_map = self.detmap.copy().astype(float)
        in_map -= bg_clip[1]
        #seg_thresh = skimage.filters.threshold_otsu(in_map, nbins=1024)
        seg_thresh = 2.0 * bg_clip[0]

        # Use photutils to segment map
        seg_map = photutils.detect_sources(in_map, threshold=seg_thresh, npixels=area_thresh, connectivity=8).array
        seg_map = AstroCell.Process.LabelShuffle(seg_map)

        # Record attributes
        self.thresh_segmap = seg_map
        self.thresh_level = seg_thresh
        self.thresh_area = area_thresh



    def WaterBorders(self, iter_total=500, seg_map=None):
        """ A method that uses a Monte Carlo series of watershed segmentations to deblend segmented cell features """

        # If no segment map specified, use map from thesholding segmentation
        if seg_map==None:
            seg_map = self.thresh_segmap

         # If segmentation map contains no segments, return null results
        if np.max(seg_map) == 0:
            self.thresh_segmap = seg_map
            return

        # Prepare parameters for Monte Carlo segmenations
        self.water_iter = iter_total
        processes = mp.cpu_count()-1

        """# Filter detection map with a Mexican-hat kernel
        kernel = astropy.convolution.kernels.Tophat2DKernel(2.0)
        self.hatmap = astropy.convolution.convolve_fft(self.detmap, kernel, interpolate_nan=True, boundary='reflect')"""

        # Run random iterations in parallel, for speed
        water_map_list = []
        pool = mp.Pool(processes=processes)
        for i in range(0, iter_total):
            if self.parallel:
                water_map_list.append( pool.apply_async( AstroCell.Process.WaterWrapper, args=(self, seg_map, iter_total,) ) )
            else:
                water_map_list.append( AstroCell.Process.WaterWrapper(copy.deepcopy(self), seg_map, iter_total) )
        pool.close()
        pool.join()
        water_map_list = [output.get() for output in water_map_list]

        # Convert walker maps to boundry maps
        border_map = np.zeros(seg_map.shape)
        for i in range(0, len(water_map_list)):
            border_map += skimage.segmentation.find_boundaries(water_map_list[i], connectivity=2)

        # Neaten border edges
        border_map[ np.where(self.thresh_segmap==0) ] = 0

        # Record watershed output
        del(water_map_list)
        gc.collect()
        self.water_border = border_map
        #astropy.io.fits.writeto('/home/chris/border_map.fits', border_map, clobber=True)



    def WalkerBorders(self, iter_total=250, seg_map=None):
        """ A method that uses a Monte Carlo series of random water segmentations to deblend segmented cell features """

        # If no segment map specified, use map from thesholding segmentation
        if seg_map==None:
            seg_map = self.thresh_segmap

         # If segmentation map contains no segments, return null results
        if np.max(seg_map) == 0:
            self.thresh_segmap = seg_map
            return

        # Prepare parameters for Monte Carlo segmenations
        self.walker_iter = iter_total
        processes = int(0.5*mp.cpu_count())

        # Run random iterations in parallel, for speed
        walker_map_list = []
        pool = mp.Pool(processes=processes)
        for i in range(0, iter_total):
            if self.parallel:
                walker_map_list.append( pool.apply_async( AstroCell.Process.WalkerWrapper, args=(self, seg_map, iter_total,) ) )
            else:
                walker_map_list.append( AstroCell.Process.WalkerWrapper(copy.deepcopy(self), seg_map, iter_total) )
        pool.close()
        pool.join()
        walker_map_list = [output.get() for output in walker_map_list]

        # Convert walker maps to boundry maps
        border_map = np.zeros(seg_map.shape)
        for i in range(0, len(walker_map_list)):
            border_map += skimage.segmentation.find_boundaries(walker_map_list[i], connectivity=2)

        # Neaten border edges
        border_map[ np.where(self.thresh_segmap==0) ] = 0

        # Record watershed output
        del(walker_map_list)
        gc.collect()
        self.walker_border = border_map
        #astropy.io.fits.writeto('/home/chris/border_map.fits', border_map, clobber=True)



    def DeblendSegment(self):
        """ Method that performs segmentation using output of watershed segmentations """

        # Perform hysteresis thresholding on this channel's watershed border map
        hyster_border = AstroCell.Process.HysterThresh(self.water_border.copy(), (0.2*self.water_iter), (0.4*self.water_iter))

        # Perform segmentation using hysteresis
        hyster_seg_map = np.invert(hyster_border).astype(int)
        hyster_seg_map[ np.where(self.thresh_segmap==0) ] = 0
        hyster_seg_map = scipy.ndimage.measurements.label(hyster_seg_map)[0]

        # Conduct binary opening to remove any remaining 'bridges'
        open_structure = scipy.ndimage.generate_binary_structure(2,1)
        hyster_seg_map_open = scipy.ndimage.binary_opening(hyster_seg_map, structure=open_structure, iterations=1).astype(float)
        hyster_seg_map_open *= hyster_seg_map

        # Remove spuriously small features, based on having 5 or fewer pixels
        hyster_seg_areas = np.unique(hyster_seg_map, return_counts=True)[1]
        hyster_exclude = np.arange(0,hyster_seg_areas.size)[ np.where(hyster_seg_areas<=5) ]
        hyster_seg_flat = hyster_seg_map.copy().flatten()
        hyster_seg_flat[np.in1d(hyster_seg_flat,hyster_exclude)] = 0
        hyster_seg_map = np.reshape(hyster_seg_flat, hyster_seg_map.shape)

        # Shuffle labels of hysteresis segmentation map, and record
        hyster_seg_map = scipy.ndimage.measurements.label(hyster_seg_map)[0]
        hyster_seg_map = AstroCell.Process.LabelShuffle(hyster_seg_map).astype(float)
        self.hyster_segmap = hyster_seg_map.copy()
        #astropy.io.fits.writeto('/home/chris/hyster_seg_map.fits', hyster_seg_map.astype(float), clobber=True)

