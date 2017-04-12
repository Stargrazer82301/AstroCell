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
in_dir = '/home/chris/Data/AstroCell/Flourescant/Liver/APCFLOX1668'#'/home/chris/Data/AstroCell/Histochemial/3100_zeb1/'#'/home/chris/Data/AstroCell/Histochemial/3100_zeb1/'#'/home/chris/Data/AstroCell/Flourescant/Mammary/Ref_LO'
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



    sdfdsfdsdsfdsvds

    astropy.io.fits.writeto('/home/chris/det_map.fits', rgb.b.detmap, clobber=True)
    astropy.io.fits.writeto('/home/chris/thresh_seg_map.fits', rgb.b.thresh_segmap, clobber=True)



    # Decide how many markers to generate, based on number of Canny and threshold segmentation features found
    n_markers = int( 5.0 * np.max([ np.unique(rgb.b.canny_features).shape[0], np.unique(rgb.b.thresh_segmap).shape[0] ]) )

    # Generate marker coordinates
    markers = np.random.random(size=(n_markers,2))
    markers[:,0] *= rgb.b.map.shape[0]
    markers[:,1] *= rgb.b.map.shape[1]

    # Prune markers located too close to each other
    markers = AstroCell.Process.ProximatePrune(markers, 0.5*np.sqrt(rgb.b.thresh_area/np.pi))











# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))
#astropy.io.fits.writeto('/home/chris/coadd.fits', coadd, clobber=True)

# Jubiliate
print('All done!')
