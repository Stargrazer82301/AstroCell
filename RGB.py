# Import smorgasbord
import os
import pdb
import dill
import multiprocessing as mp
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
import joblib
import astrodendro
import AstroCell.Process
import AstroCell.Image
from ChrisFuncs import SigmaClip
plt.ioff()




class RGB():
    """ Class that stores an RGB triplet of image arrays; the methods of an RGB object concern handling all three channels together"""



    def __init__(self, in_path):
        """ Initialise the RGB object """

        # Store path for future use
        self.path = in_path
        self.dir_name = os.path.split(in_path)[0]
        self.file_name = os.path.split(in_path)[1]

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
        dill_file = '.'.join(self.file_name.split('.')[:-1])+'.dj'
        dill_path = os.path.join(dill_dir, dill_file)

        # Conduct pickling
        dill.dump( self, open( dill_path, 'wb' ) )



    def TempDir(self, temp):
        """ Brief method that records loation of current temporary directory to the RGB object and it's constituent Image objects """

        # Record temp dir location attribute to self
        self.temp = temp

        # Record temp dir location attribute to Image objects
        for channel in self.iter:
            channel.temp = temp



    def RecParallel(self, parallel):
        """ Brief method that records whether utilising multiprocessing to the RGB object and it's constituent Image objects """

        # Record temp dir location attribute to self
        self.parallel = parallel

        # Record temp dir location attribute to Image objects
        for channel in self.iter:
            channel.parallel = parallel



    def MakeCoadd(self):
        """ Method that adds an new Image object, that's a coadd of the three channels """

        # Initialise AstroCell.Image object for a coadd of all three channels
        coadd = np.sum(self.cube,axis=2).astype(float) / 3.0
        self.coadd = AstroCell.Image.Image(coadd)

        # Create tuple containing each channel's Image object (including the coadd) for iterating over
        self.iter_coadd = (self.b,self.g,self.r,self.coadd)

        # Make tuple of strings giving name of each channel (including the coadd), to use for iteration
        self.channels = ('b','g','r','coadd')

        # Record name of coadd channel, location of temp dir, and parallelisation
        self.coadd.name = 'coadd'
        self.coadd.temp = self.temp
        self.coadd.parallel = self.parallel



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
        if self.parallel:
            logdog_blob_list = joblib.Parallel( n_jobs=mp.cpu_count()-1 )\
                                              ( joblib.delayed( channel.LogDogBlobs )\
                                              ( canny_features=self.coadd.canny_features )\
                                              for channel in self.iter_coadd )
            [ setattr(self.iter_coadd[c],'logdog_mask',logdog_blob_list[c][0]) for c in range(0,len(self.iter_coadd)) ]
            [ setattr(self.iter_coadd[c],'logdog_features',logdog_blob_list[c][1]) for c in range(0,len(self.iter_coadd)) ]

        # Else if not operating in parallel, do things the straightforward way
        else:
            [ channel.LogDogBlobs(canny_features=self.coadd.canny_features) for channel in self.iter_coadd ]



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



    def DetFilter(self):
        """ Method that removes smooth all large-scale background structure from map, to create an optimal detection map """

        # Use Canny features map to get distribution of feature sizes (assuming circuar features)
        canny_areas = np.unique(self.coadd.canny_features, return_counts=True)[1].astype(float)
        canny_diams = 2.0 * np.sqrt( canny_areas / np.pi)
        #astropy.io.fits.writeto('/home/chris/canny_coadd.fits', canny_coadd, clobber=True)

        # Decide size of filter to apply, based upon typical size range of Canny cells
        canny_diams_clip = SigmaClip(canny_diams, median=True)
        kernel_size = np.percentile(canny_diams_clip[1]+canny_diams_clip[0], 90.0)
        kernel = astropy.convolution.kernels.Gaussian2DKernel(kernel_size)

        # Iterate over each channel, creating and removing background model for each
        for channel in self.iter_coadd:

            # Create background map by excluding Canny cells from image
            canny_bg_map = channel.map.copy().astype(float)
            canny_bg_map[ np.where(self.blob_mask>0) ] = np.NaN
            canny_bg_fill = SigmaClip(canny_bg_map, median=True, sigma_thresh=3.0)[1]

            # Also exclude aberantly bright pixels from background map
            canny_bg_clip = SigmaClip(canny_bg_map, median=True, sigma_thresh=5.0)
            canny_bg_map[ np.where(canny_bg_map>(canny_bg_clip[1]+(5.0*canny_bg_clip[0]))) ] = np.NaN

            # Smooth background mapusing a Gaussian filter
            conv_map = astropy.convolution.convolve_fft(canny_bg_map, kernel, interpolate_nan=True, normalize_kernel=True, quiet=True,
                                                        boundary='reflect', fill_value=canny_bg_fill, allow_huge=True)
            conv_map[ np.where( np.isnan(channel.map)==True ) ] = np.NaN

            # Subtract smoothed background map from original image to make detmap
            conv_sub = channel.map - conv_map

            # Re-set zero level, and record map to object
            conv_sub -= np.nanmin(conv_sub)
            channel.detmap = conv_sub
            #astropy.io.fits.writeto('/home/chris/conv.fits', conv_map, clobber=True)



    def OverviewImage(self):

        """
        rand_cmap = photutils.utils.random_cmap(np.nanmax(deblend_map.array) + 1, random_state=12345)
        norm = astropy.visualization.mpl_normalize.ImageNormalize( stretch=astropy.visualization.AsinhStretch() )
        map_aspect = float(rgb.coadd.map.shape[1]) / float(rgb.coadd.map.shape[0])
        fig_x_panes = 3
        fig_y_panes = 1
        fig_aspect = ( fig_x_panes * map_aspect ) / fig_y_panes
        fig_x_dim = 10.0 * fig_aspect
        fig_y_dim = 10.0
        fig, axes = plt.subplots(fig_y_panes, fig_x_panes, figsize=(fig_x_dim, fig_y_dim))
        axes[0].imshow(in_map, origin='lower', cmap='inferno', vmin=np.percentile(in_map,5), vmax=np.percentile(in_map, 95))
        axes[1].imshow(seg_map, origin='lower', cmap=rand_cmap)
        axes[2].imshow(deblend_map.array, origin='lower', cmap=rand_cmap)
        [ ax.set_xticklabels([]) for ax in axes ]
        [ ax.set_yticklabels([]) for ax in axes ]
        fig.tight_layout()
        fig.savefig( os.path.join( out_dir, in_image.replace('.bmp','_phoutils_seg.png') ), dpi=400.0 )
        """



class Raw(object):
    """ A dummy class, to allow 'holding' objects to be created """
    pass

