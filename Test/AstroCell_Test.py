# Import smorgasbord
import pdb
import sys
sys.path.append('/home/chris/Dropbox/Work/Scripts/')
import os
import shutil
import imghdr
import inspect
import time
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import scipy.stats
import scipy.ndimage
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
import skimage.restoration
import PIL.Image
import AstroCell
import AstroCell.RGB
import AstroCell.Image
import AstroCell.IO
plt.ioff()

# Include reloads, to handle changes
import importlib
importlib.reload(AstroCell)
importlib.reload(AstroCell.RGB)
importlib.reload(AstroCell.Image)
importlib.reload(AstroCell.IO)



# State input directory and create output directory inside it
in_dir = '/home/chris/Data/AstroCell/Histochemial/3100_zeb1/'#'/home/chris/Data/AstroCell/Flourescant/Liver/APCFLOX1668'#'/home/chris/Data/AstroCell/Histochemial/3100_zeb1/'#'/home/chris/Data/AstroCell/Flourescant/Mammary/Ref_LO'
out_dir = os.path.join(in_dir, 'AstroCell_Output')
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

# Initialise temp directory class
temp = AstroCell.IO.TempDir(out_dir)

# Identify and loop over all image files in input directory
in_files = os.listdir(in_dir)
in_files = [in_file for in_file in in_files if not os.path.isdir(os.path.join(in_dir,in_file))]
in_images = [in_file for in_file in in_files if imghdr.what(os.path.join(in_dir,in_file))!=None]
for in_image in np.random.permutation(in_images):



    # Read in raw image, constructing an AstroCell RGB object
    rgb = AstroCell.RGB.RGB(os.path.join(in_dir, in_image))

    # Clean edges of images
    [ channel.CleanEdges() for channel in rgb.iter ]

    # Create coadd of all three channels
    rgb.MakeCoadd()

    # Craete Canny-based mask identifying pixels that contain cells
    rgb.CannyMask()

    # Determine if image is black-background; if not, set it so that it is
    rgb.BlackOnWhite()

    # Remove large-scale background structures from image (to create source extraction map)
    rgb.DetFilter()

    """# Construct basic matched filter in each channel, using Canny features
    [ channel.CannyCellStack() for channel in rgb.iter_coadd ]"""

    # Use canny features to create markers for cells and background, to anchor segmentation
    [ channel.ThreshSegment(rgb.canny_mask) for channel in rgb.iter_coadd ]




    """
    from ChrisFuncs import SigmaClip
    in_map = rgb.r.detmap.copy().astype(float)
    bg_map = in_map.copy().astype(float)
    bg_map[ np.where(rgb.r.canny_features>0) ] = np.NaN
    clip = SigmaClip(bg_map, median=True, sigma_thresh=5.0)
    in_map -= clip[1]
    thresh = 1.0 * clip[0]

    #thresh = skimage.filters.threshold_otsu(in_map, nbins=1024)

    canny_areas = np.unique(rgb.r.canny_features, return_counts=True)[1].astype(float)
    area_thresh = int(np.round(np.percentile(canny_areas, 10.0)))

    seg_map = photutils.detect_sources(in_map, threshold=thresh, npixels=area_thresh, connectivity=4)
    seg_areas = np.unique(seg_map, return_counts=True)[1].astype(float)

    conv_map = astropy.convolution.convolve_fft(in_map, rgb.r.canny_stack, interpolate_nan=True, normalize_kernel=True, boundary='reflect', allow_huge=True)




    tophat_kernel = astropy.convolution.Tophat2DKernel( (np.median(canny_areas)/np.pi)**0.5 )
    #deblend_map = photutils.deblend_sources(in_map, seg_map, npixels=area_thresh, filter_kernel=rgb.r.canny_stack, labels=None, nlevels=128, contrast=0.0001, mode='exponential', connectivity=4, relabel=True)
    deblend_map = photutils.deblend_sources(in_map, seg_map, npixels=area_thresh, filter_kernel=tophat_kernel, labels=None, nlevels=1024, contrast=0.000000001, mode='exponential', connectivity=8, relabel=True)

    grey_erode_structure = scipy.ndimage.generate_binary_structure(2,1)
    grey_erode_map = scipy.ndimage.grey_erosion(in_map, structure=grey_erode_structure, mode='reflect')
    astropy.io.fits.writeto('/home/chris/red_grey_erode.fits', grey_erode_map, clobber=True)


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




    sdfdsfds
    pdb.set_trace()











# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))
#astropy.io.fits.writeto('/home/chris/coadd.fits', coadd, clobber=True)

# Jubiliate
print('All done!')
