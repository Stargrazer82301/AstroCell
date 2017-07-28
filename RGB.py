# Import smorgasbord
import os
import pdb
import shutil
import dill
import warnings
warnings.simplefilter('ignore', category=Warning)
import multiprocessing as mp
import numpy as np
import scipy.stats
import matplotlib.pylab as plt
import matplotlib.image
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.convolution
import astropy.stats
import astropy.visualization
import astropy.visualization.mpl_normalize
import astropy.io.fits
import astropy.table
import photutils
import skimage.feature
import sklearn.cluster
import PIL.Image
import joblib
import webcolors
import AstroCell.Process
import AstroCell.Image
from ChrisFuncs import SigmaClip
plt.ioff()




class RGB():
    """ Class that stores an RGB triplet of image arrays; the methods of an RGB object concern handling all three channels together"""



    def __init__(self, in_path, out_dir):

        """ Initialise the RGB object """
        # Store path for future use
        self.path = in_path
        self.in_dir = os.path.split(in_path)[0]
        self.out_dir = out_dir
        self.in_file = os.path.split(in_path)[1]

        # Report processing started for this image
        self.id = self.in_file.replace(' ','_')

        # Read in the image file, and convert to array
        bitmap_image = PIL.Image.open(in_path)
        rgb_image = np.array(bitmap_image)

        # Store the full rgb cube, as a float
        self.cube = rgb_image.astype(float)

        # Initialise primary AstroCell.Image object for each channel
        self.r = AstroCell.Image.Image(self.cube[:,:,0].copy())
        self.g = AstroCell.Image.Image(self.cube[:,:,1].copy())
        self.b = AstroCell.Image.Image(self.cube[:,:,2].copy())

        # Create tuple containing each channel's Image object for iterating over
        self.iter = (self.b,self.g,self.r)

        # Make tuple of strings giving name of each channel, to use for iteration
        self.channels = ('b','g','r')

        # Set data type for each channel to be float
        for channel in self.iter:
            channel.map = channel.map.astype(float)

        # Record each channel's name as an attribute
        for i in range(0, len(self.iter)):
            self.iter[i].name =  self.channels[i]



    def Dill(self, dill_dir):
        """ A method that dill pickles an RGB object, so that it can re-used later, to save reprocessing (especially during testing) """

        # Construct dill pickle's file name
        dill_file = '.'.join(self.in_file.split('.')[:-1])+'.dj'
        dill_path = os.path.join(dill_dir, dill_file)

        # Conduct pickling
        dill.dump( self, open( dill_path, 'wb' ) )



    def TempDir(self):
        """ Brief method that creates a temporary directory, and records the path to the RGB object and constituent Image objects """

        # Record temp dir location attribute to self
        self.temp_dir = os.path.join(self.out_dir,'Temp')
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
        os.mkdir(self.temp_dir)

        # Record temp dir location attribute to Image objects
        for channel in self.iter:
            channel.temp_dir = self.temp_dir



    def TempDirTidy(self):
        """ Method that deletes temporary directory at end of processing """

        # Record temp dir location attribute to self
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)



    def RecVerbose(self, verbose):
        """ Brief method that records verbosity each Image object """

        # Record parallel status attribute to self
        self.verbose = verbose

        # Record parallel status attribute to Image objects
        for channel in self.iter:
            channel.verbose = self.verbose

        # If verbosity requested, report that processing has commenced on current image file
        if self.verbose:
            print('['+self.id+'] Loading image '+self.in_file+'.')



    def RecParallel(self, parallel):
        """ Brief method that records whether utilising multiprocessing to the RGB object and it's constituent Image objects """

        # Record parallel status attribute to self
        self.parallel = parallel

        # Record parallel status attribute to Image objects
        for channel in self.iter:
            channel.parallel = self.parallel



    def RecMCFactor(self, mc_factor):
        """ Brief method that records what multiplier is to be used to scale number of iterations for Monte-Carlo processes """

        # Record MC factor attribute to self
        self.mc_factor = mc_factor

        # Record MC factor attribute to Image objects
        for channel in self.iter:
            channel.mc_factor = mc_factor



    def MakeCoadd(self):
        """ Method that adds an new Image object, that's a coadd of the three channels """

        # Initialise AstroCell.Image object for a coadd of all three channels
        coadd = np.sum(self.cube,axis=2).astype(float) / 3.0
        self.coadd = AstroCell.Image.Image(coadd)

        # Create tuple containing each channel's Image object (including the coadd) for iterating over
        self.iter_coadd = (self.coadd,self.b,self.g,self.r)

        # Make tuple of strings giving name of each channel (including the coadd), to use for iteration
        self.channels = ('coadd','b','g','r')

        # Record name of coadd channel, location of temp dir, and parallelisation
        self.coadd.name = 'coadd'
        self.coadd.verbose = self.verbose
        self.coadd.parallel = self.parallel
        self.coadd.mc_factor = self.mc_factor
        self.coadd.temp_dir = self.temp_dir



    def CannyMask(self):
        """ A method that combines the Canny blob detectionn in each channel to create a mask tracing cell-filled regions """

        # Do a rough Canny-based cell extraction on each channel, including the coadd
        self.canny_cube = np.zeros([self.cube.shape[0],self.cube.shape[1],4])
        for i in range(0, len(self.iter_coadd)):
            self.canny_cube[:,:,i] = self.iter_coadd[i].canny_features

        # Coadd the Canny feature maps from each channel
        canny_coadd = np.sum(self.canny_cube, axis=2)
        canny_where = np.where(canny_coadd>0)
        #astropy.io.fits.writeto('/home/chris/canny_coadd.fits', canny_coadd, clobber=True)

        # Create Canny mask, and record
        canny_mask = canny_coadd.copy()
        canny_coadd[canny_where] = 1
        self.canny_mask = canny_mask



    def BlackOnWhite(self):
        """ A method that determines whether an image is on a black background; if not, image is inverted so that it is """

        # Flag pixel values for pixels within Canny features
        cube = self.cube.copy().astype(int)
        canny_where = np.where(self.canny_mask>0)
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



    def LogDogBlobsWrapper(self):
        """ Wrapper around the optionally-parallel LoG-DoG blob finding Image method """

        # If parallel operation requested, process channels simultaneously using joblib
        if self.parallel.parallel:
            logdog_blob_list = joblib.Parallel( n_jobs=self.parallel.subthreads )\
                                              ( joblib.delayed( channel.LogDogBlobs )\
                                              ( canny_features=self.coadd.canny_features )\
                                              for channel in self.iter_coadd )
            [ setattr(self.iter_coadd[c],'logdog_mask',logdog_blob_list[c][0]) for c in range(0,len(self.iter_coadd)) ]
            [ setattr(self.iter_coadd[c],'logdog_features',logdog_blob_list[c][1]) for c in range(0,len(self.iter_coadd)) ]

        # Else if not operating in parallel, do things the straightforward way
        else:
            [ channel.LogDogBlobs(canny_features=self.coadd.canny_features, force_attribute=True) for channel in self.iter_coadd ]



    def BlobMask(self):
        """ A method that combines the Canny, LoG, and DoG blob detectionn in each channel to create a mask of cell-filled regions """

        # Create cube of Canny features from each channel (including coadd)
        self.canny_cube = np.zeros([self.cube.shape[0],self.cube.shape[1],4])
        for i in range(0, len(self.iter_coadd)):
            self.canny_cube[:,:,i] = self.iter_coadd[i].canny_features

        # Create cube of LoG-DoG blobs from each channel (including coadd)
        logdog_cube = np.zeros([self.cube.shape[0],self.cube.shape[1],4])
        for i in range(0, len(self.iter_coadd)):
            logdog_cube[:,:,i] = self.iter_coadd[i].logdog_features

        # Coadd the Canny feature maps and LoG-DoG blob from each channel
        canny_coadd = np.sum(self.canny_cube, axis=2)
        canny_coadd[ np.where( canny_coadd == np.nanmin(canny_coadd) ) ] = 0
        logdog_coadd = np.sum(logdog_cube, axis=2)
        blob_coadd = canny_coadd + logdog_coadd
        #astropy.io.fits.writeto('/home/chris/blob_coadd.fits', blob_coadd, clobber=True)

        # Create Canny mask, and record
        blob_mask = blob_coadd.copy()
        blob_mask[np.where(blob_coadd>0)] = 1
        blob_mask[np.where(blob_coadd<=0)] = 0
        self.blob_mask = blob_mask

        # Create per-channel blob masks
        for i in range(0, len(self.iter_coadd)):
            channel_blob_coadd = self.canny_cube[:,:,i] + logdog_cube[:,:,i]
            channel_blob_mask = channel_blob_coadd.copy()
            channel_blob_mask[np.where(channel_blob_coadd>0)] = 1
            channel_blob_mask[np.where(channel_blob_coadd<=0)] = 0
            self.iter_coadd[i].blob_mask = channel_blob_mask



    def DetFilter(self):
        """ Method that removes smooth all large-scale background structure from map, to create an optimal detection map """

        # Use Canny features map to get distribution of feature sizes (assuming circuar features)
        canny_areas = np.unique(self.coadd.canny_features, return_counts=True)[1].astype(float)
        canny_diams = 2.0 * np.sqrt( canny_areas / np.pi)
        #astropy.io.fits.writeto('/home/chris/canny_coadd.fits', canny_coadd, clobber=True)

        # Decide size of filter to apply, based upon typical size range of Canny cells
        canny_diams_clip = SigmaClip(canny_diams, median=True)
        kernel_size = 0.5 * canny_diams_clip[1] #np.percentile(canny_diams_clip[1]+canny_diams_clip[0], 90.0)
        kernel = astropy.convolution.kernels.Gaussian2DKernel(kernel_size)

        # Iterate over each channel, creating and removing background model for each
        for channel in self.iter_coadd:

            # Create background map by excluding Canny cells from image
            canny_bg_map = channel.map.copy().astype(float)
            canny_bg_map[ np.where(channel.blob_mask>0) ] = np.NaN
            canny_bg_fill = SigmaClip(canny_bg_map, median=True, sigma_thresh=3.0)[1]

            # Also exclude aberantly bright pixels from background map
            canny_bg_clip = SigmaClip(canny_bg_map, median=True, sigma_thresh=5.0)
            canny_bg_map[ np.where(canny_bg_map>(canny_bg_clip[1]+(5.0*canny_bg_clip[0]))) ] = np.NaN

            # Smooth background mapusing a Gaussian filter
            conv_map = astropy.convolution.convolve_fft(canny_bg_map, kernel, nan_treatment='interpolate', normalize_kernel=True, quiet=True,
                                                        boundary='wrap', fill_value=canny_bg_fill, allow_huge=True)
            conv_map[ np.where( np.isnan(channel.map)==True ) ] = np.NaN

            # Subtract smoothed background map from original image to make detmap
            conv_sub = channel.map - conv_map

            # Re-set zero level, and record map to object
            conv_sub -= np.nanmin(conv_sub)
            channel.bgmap = conv_map
            channel.detmap = conv_sub
            #astropy.io.fits.writeto('/home/chris/conv.fits', conv_map, clobber=True)



    def SegmentCombine(self):
        """ Method that combines the segmentations from each individual channel to produces the final segmenation """

        # Create segmentation 'cube', holding the segments from each band (NB, the channel index comes first, to simplify FITS output)
        hyster_seg_cube = np.zeros([4, self.cube.shape[0], self.cube.shape[1]]).astype(int)
        for i in range(0, len(self.iter_coadd)):
            hyster_seg_cube[i,:,:] = self.iter_coadd[i].hyster_segmap

        # Create segmentation stack
        hyster_seg_stack = np.sum(hyster_seg_cube.astype(bool).astype(int), axis=0)

        # Work out combined minimum area threshold
        thresh_area_list = []
        [ thresh_area_list.append(channel.thresh_area) for channel in self.iter_coadd ]
        self.thresh_area = self.coadd.thresh_area#np.nanmin(np.array(thresh_area_list))

        # Convolve segmentation stack with a small Gaussian kernel (but keep empty pixels as empty)
        hyster_seg_stack_mask = np.zeros(hyster_seg_stack.shape).astype(int)
        hyster_seg_stack_mask[np.where(hyster_seg_stack>0)] = 1
        kernel = astropy.convolution.kernels.Tophat2DKernel(1.5)
        hyster_seg_stack = astropy.convolution.convolve_fft(hyster_seg_stack, kernel, nan_treatment='interpolate', normalize_kernel=True,
                                                            quiet=True, boundary='reflect', fill_value=0, allow_huge=True)
        hyster_seg_stack[np.where(hyster_seg_stack_mask==0)] = 0

        # Create stacked thresholding segmentation map from all 4 channels, to identify every pixel deemed to be 'cell'
        thresh_seg_cube = np.zeros([4, self.cube.shape[0], self.cube.shape[1]]).astype(int)
        for i in range(0, len(self.iter_coadd)):
            thresh_seg_cube[i,:,:] = self.iter_coadd[i].thresh_segmap
        thresh_seg_stack = np.sum(thresh_seg_cube.astype(bool).astype(int), axis=0)

        # Initiate meta-segmentation as its own (quasi-dummy) Image object
        self.meta = AstroCell.Image.Image(hyster_seg_stack)
        self.meta.name = 'meta'
        self.meta.verbose = self.verbose
        self.meta.parallel = self.parallel
        self.meta.mc_factor = self.mc_factor
        self.meta.temp_dir = self.temp_dir
        self.meta.detmap = self.meta.map.copy()
        self.meta.thresh_area = self.thresh_area
        self.meta.thresh_segmap = thresh_seg_stack
        #self.meta.thresh_segmap = scipy.ndimage.measurements.label(hyster_seg_stack.astype(bool).astype(int))[0]

        # Do blob-finding required for later segmentation
        self.meta.CannyBlobs(sigma=2.0)
        self.meta.LogDogBlobs(canny_features=None, force_attribute=True)

        # Perform meta-segmentation using monte-carlo watershed thresholding
        self.meta.WaterBorders(seg_map=self.meta.thresh_segmap)
        self.meta.water_border[np.where(self.meta.thresh_segmap==0)] = 0

        # Perform hysteresis thresholding on watershed border map
        self.meta.DeblendSegment(thresh_lower=0.4, thresh_upper=0.9, meta=True)

        # Record final segmentation map to object
        self.segmap = self.meta.hyster_segmap.astype(int)



    def CellPhotom(self):
        """ Method that measures the properties of all the segmented cells in the image """

        # Declare columns for 'photometry' table
        table_col_names = ('id',
                           'i_coord','j_coord','area',
                           'b_flux','g_flux','r_flux',
                           'b_mu','g_mu','r_mu',
                           'bg_colour','br_colour','gr_colour',
                           'b_top_mu','g_top_mu','r_top_mu',
                           'bg_top_colour','br_top_colour','gr_top_colour')
        table_col_dtypes = ('i8',
                            'f8','f8','f8',
                            'f8','f8','f8',
                            'f8','f8','f8',
                            'f8','f8','f8',
                            'f8','f8','f8',
                            'f8','f8','f8')

        # Create table object to hold cell properties
        self.table = astropy.table.Table(names=table_col_names, dtype=table_col_dtypes)

        # Identify some specific sub-sets of the columns
        self.table.data_cols = ('b_flux','g_flux','r_flux',
                                'b_mu','g_mu','r_mu',
                                'bg_colour','br_colour','gr_colour',
                                'b_top_mu','g_top_mu','r_top_mu',
                                'bg_top_colour','br_top_colour','gr_top_colour')
        self.table.colour_cols = ('bg_top_colour','br_top_colour','gr_top_colour')

        # Loop over segments
        for i in range(1, int(self.segmap.max())):

            # Skip non-existant segments
            seg_where = np.where(self.segmap==i)
            if seg_where[0].shape[0] == 0:
                continue

            # Calculate area and coordinates
            row_area = seg_where[0].shape[0]
            index_map = self.segmap.copy()
            index_map[np.where(self.segmap!=i)] == 0
            row_coord_i = np.median( np.where(self.segmap==i)[0] )
            row_coord_j = np.median( np.where(self.segmap==i)[0] )

            # Calculate fluxes
            row_flux_b = np.sum(self.b.map[seg_where])
            row_flux_g = np.sum(self.g.map[seg_where])
            row_flux_r = np.sum(self.r.map[seg_where])

            # Calculate mean surface brightnesses
            row_mu_b = np.mean(self.b.map[seg_where])
            row_mu_g = np.mean(self.g.map[seg_where])
            row_mu_r = np.mean(self.r.map[seg_where])

            # Calculate mean-mu colours (with offset of 10 to prevent zero-based log troubles)
            row_colour_bg = np.log10(10.0+row_mu_b) - np.log10(10.0+row_mu_g)
            row_colour_br = np.log10(10.0+row_mu_b) - np.log10(10.0+row_mu_r)
            row_colour_gr = np.log10(10.0+row_mu_g) - np.log10(10.0+row_mu_r)

            # Calculate top-mu colours
            med_mu_b = np.median(self.b.map[seg_where])
            med_mu_g = np.median(self.g.map[seg_where])
            med_mu_r = np.median(self.r.map[seg_where])
            top_where_b = np.where( (self.segmap==i) & (self.b.map>med_mu_b) )
            top_where_g = np.where( (self.segmap==i) & (self.g.map>med_mu_g) )
            top_where_r = np.where( (self.segmap==i) & (self.r.map>med_mu_r) )
            row_top_mu_b = np.nanmax([ 0.0, np.mean(self.b.map[top_where_b]) ])
            row_top_mu_g = np.nanmax([ 0.0, np.mean(self.g.map[top_where_g]) ])
            row_top_mu_r = np.nanmax([ 0.0, np.mean(self.r.map[top_where_r]) ])
            row_top_colour_bg = np.log10(10.0+row_top_mu_b) - np.log10(10.0+row_top_mu_g)
            row_top_colour_br = np.log10(10.0+row_top_mu_b) - np.log10(10.0+row_top_mu_r)
            row_top_colour_gr = np.log10(10.0+row_top_mu_g) - np.log10(10.0+row_top_mu_r)

            # Calculate and record segment properties
            self.table.add_row([i,
                                row_area, row_coord_i, row_coord_j,
                                row_flux_b, row_flux_g, row_flux_r,
                                row_mu_b, row_mu_g, row_mu_r,
                                row_colour_bg, row_colour_br, row_colour_gr,
                                row_top_mu_b,  row_top_mu_g,  row_top_mu_r,
                                row_top_colour_bg, row_top_colour_br,row_top_colour_gr])



    def CellClassify(self, cell_colours=None):
        """ Method that classifies cells based on their location in a multi-dimensional parameter space """

        # Select only table columns that contain data, convert them into an array for ease of processing
        data_array = self.table[self.table.colour_cols].as_array()
        data_array = data_array.view((float, len(data_array.dtype.names)))

        # Scale data for clustering analysis
        data_scaler = sklearn.preprocessing.StandardScaler(copy=False)
        data_array = data_scaler.fit_transform(data_array)

        # If number os cell colours is pre-stated, use Ward Agglomerative custer-finding algorithm
        if isinstance(cell_colours, (float,int)):
            cluster_algorithm = sklearn.cluster.AgglomerativeClustering(n_clusters=cell_colours)
            cluster_algorithm.fit(data_array)

        # If number of cell colours not pre-stated, then use mean-shift cluster-finding algorithm
        elif cell_colours == None:
            cluster_bandwidth = sklearn.cluster.estimate_bandwidth(data_array)
            cluster_algorithm = sklearn.cluster.MeanShift(bandwidth=cluster_bandwidth, n_jobs=-2)

        # Fit cluster-finding algorthm to data, and record claifications to table and object
        cluster_algorithm.fit(data_array)
        self.table['label'] = cluster_algorithm.labels_
        self.labels_list =  self.table['label'].tolist()
        self.labels = np.unique(self.labels_list)

        # Record label counts
        self.labels_counts = []
        for i in self.labels:
            self.labels_counts.append(np.array(self.labels_list)[ np.where( np.array(self.labels_list) == self.labels[i] ) ].size)
        self.labels_counts = np.array(self.labels_counts).astype(int)

        # Create classification map of image
        self.labelmap = np.zeros(self.coadd.map.shape) - 1
        for i in range(0,len(self.table)):
            self.labelmap[ np.where( self.segmap == int(self.table['id'][i]) ) ] = self.table['label'][i]

        """# Create label maps of image, dilated slightly for display purposes
        self.labelcube = np.zeros([ np.unique(self.table['label']).size, self.coadd.map.shape[0], self.coadd.map.shape[1] ])"""

        # Write 'photometry' table to output directory
        self.table.write(os.path.join(self.out_dir,'.'.join(self.id.split('.')[:-1])+'.csv'))



    def CellColours(self):
        """ Method that works out the colours of classified cells """

        # Take reference colour names from webcolors package, and convert to hsv to obtain hue
        wc_colours_hex = webcolors.css3_names_to_hex
        wc_colours_names = np.array(list(webcolors.css3_names_to_hex.keys()))
        wc_colours_hsv = np.zeros([ len(wc_colours_hex.keys()), 3 ])
        wc_colours_lab = np.zeros([ len(wc_colours_hex.keys()), 3 ])
        for i in range(0, len(wc_colours_names)):
            wc_colour_hex = wc_colours_hex[wc_colours_names[i]]
            wc_colour_rgb = np.array(list(webcolors.hex_to_rgb(wc_colour_hex)))
            wc_colours_hsv[i,:] = skimage.color.rgb2hsv([[wc_colour_rgb]])[0][0]
            wc_colours_lab[i,:] = skimage.color.rgb2lab([[wc_colours_hsv[i,:]]])[0][0]

        # Sort refence colours in order of hue
        wc_colours_argsort = np.argsort(wc_colours_hsv[:,0])
        wc_colours_names = wc_colours_names[wc_colours_argsort]
        wc_colours_hsv = wc_colours_hsv[wc_colours_argsort]
        wc_colours_lab = wc_colours_lab[wc_colours_argsort]

        # Create storage variables for label colours
        labels_colours_dist = np.zeros([ len(np.unique(self.labels_list)), len(wc_colours_names) ])
        labels_rgb = np.zeros([ len(np.unique(self.labels_list)), 3 ])
        labels_rgb_names = ['']*len(np.unique(self.labels_list))
        labels_hsv = np.zeros([ len(np.unique(self.labels_list)), 3 ])
        labels_rgb_bright = np.zeros([ len(np.unique(self.labels_list)), 3 ])

        # Loop over labels, in order to work out what colour they actually represent
        for i in self.labels:

            # Work out 'typical' rgb colour for pixels in cells assigned this label
            labels_rgb[i,0] = np.nanmedian( self.table['r_top_mu'] ) / 255
            labels_rgb[i,1] = np.nanmedian( self.table['g_top_mu'] ) / 225
            labels_rgb[i,2] = np.nanmedian( self.table['b_top_mu'] ) / 225

            # Convert colour from RGB to HSV (useful to generate bright 'clean' version of colour)
            label_hsv = skimage.color.rgb2hsv(np.array([[labels_rgb[i,:]]]))[0][0]
            labels_hsv[i,:] = label_hsv
            labels_rgb_bright[i,:] = skimage.color.hsv2rgb( np.array([[np.array([label_hsv[0],1,1])]]) )[0][0]

            # Now convert to LAB olour (useful, as this has actual perceptual meaning)
            label_lab = skimage.color.rgb2lab( skimage.color.hsv2rgb( np.array([[label_hsv]]) ) )[0][0]

            # Record colour distances from defined colours
            for j in range(0, len(wc_colours_names)):
                labels_colours_dist[i,j] = scipy.spatial.distance.euclidean(label_lab, wc_colours_lab[j,:])

            # Work out which webcolors reference colour is closes to this colour
            labels_rgb_names[i] = AstroCell.Process.ColourName(tuple((255*labels_rgb_bright[i,:]).tolist()))

        # Sort labels into order of hue, for consistency in output files, and record
        labels_colours_argsort = np.argsort(labels_hsv[:,0])
        self.labels = self.labels[labels_colours_argsort]
        self.labels_counts = self.labels_counts[labels_colours_argsort]
        self.labels_rgb = labels_rgb[labels_colours_argsort]
        self.labels_rgb_bright = labels_rgb_bright[labels_colours_argsort]
        self.labels_rgb_names = np.array(labels_rgb_names)[labels_colours_argsort]



    def OutFileRecord(self):
        """ Method that writes results of cell count to output text file """

        # Open the output file
        table_file_path = os.path.join(self.out_dir, 'AstroCell_Results.txt')
        table_file = open( table_file_path, 'a')

        # Construct output line string
        table_line = str(self.id)
        for i in range(0, self.labels.size):
            table_line += ','+self.labels_rgb_names[i]+','+str(self.labels_counts[i])

        # Write to file, and close
        table_line += '\n'
        table_file.write(table_line)
        table_file.close()



    def OverviewImages(self):
        """ Method that produces overview image figure of the results """

        # Use size of image, and number of panes desired, to determine size of figure
        map_aspect = float(self.coadd.map.shape[1]) / float(self.coadd.map.shape[0])
        fig_x_panes = 1
        fig_y_panes = 2
        fig_aspect = ( fig_x_panes * map_aspect ) / fig_y_panes
        fig_x_dim = 10.0 * fig_aspect
        fig_y_dim = 10.0

        # Create figure and axes
        import warnings
        warnings.filterwarnings("ignore")
        plt.ioff()
        fig, axes = plt.subplots(fig_y_panes, fig_x_panes, figsize=(fig_x_dim, fig_y_dim))

        # Remove ticks and frames from axes
        [ ax.set_xticklabels([]) for ax in axes ]
        [ ax.set_yticklabels([]) for ax in axes ]
        [ ax.axhline(linewidth=5, color='black') for ax in axes ]
        [ ax.axvline(linewidth=5, color='black') for ax in axes ]
        [ ax.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off') for ax in axes ]

        # Generate dilated verison of seg map, for display
        image_label_dilate = np.zeros(self.segmap.shape) - 1
        for i in np.unique(self.labels_list):
            for j in range(0, len(self.table)):
                if self.table['label'][j] == i:
                    image_label_indv = np.zeros(self.segmap.shape)
                    image_label_indv[ np.where(self.segmap==self.table['id'][j]) ] = 1
                    image_label_indv = scipy.ndimage.morphology.binary_erosion(image_label_indv)
                    image_label_dilate[ np.where(image_label_indv==1) ] = i
        image_label_dilate[np.where(self.segmap==0)] = -1

        # Plot label borders onto image, for each label
        image_label_border = self.cube.copy() / 255
        for i in range(0, self.labels.size):
            image_label_temp = image_label_dilate.copy()
            image_label_temp[np.where(image_label_temp!=self.labels[i])] = -1
            image_label_border = skimage.segmentation.mark_boundaries(image_label_border, image_label_temp,
                                                                      color=self.labels_rgb_bright[i,:], mode='thick', background_label=-1)

        # Shade in greyscale image regions according to label
        image_label_colour = skimage.color.label2rgb(self.labelmap, image=image_label_border,
                                                     colors=self.labels_rgb_bright[self.labels].tolist(), bg_label=-1, bg_color=[1,1,1], image_alpha=0.999)

        # Plot image panels
        axes[0].imshow(image_label_border, origin='upper')
        axes[1].imshow(image_label_colour, origin='upper')

        # Set figure to have tight layout, remove axis frames, and save to file
        fig.tight_layout()
        fig.savefig( os.path.join( self.out_dir, '.'.join(self.id.split('.')[:-1])+'_output.png' ), dpi=200 )
        #astropy.io.fits.writeto('/home/chris/bob.fits', image_label_dilate.astype(float), clobber=True)


