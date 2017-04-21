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
import multiprocessing as mp
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
import skimage.feature
import skimage.restoration
import joblib
import AstroCell
import AstroCell.RGB
import AstroCell.Image
import AstroCell.IO
plt.ioff()

# Identify location and set Dropbox path accordingly
import socket
location = socket.gethostname()
if location == 'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'sputnik':
    dropbox = '/home/chris/Dropbox/'
if location == 'saruman':
    dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Include reloads, to handle any recent changes
import importlib
importlib.reload(AstroCell)
importlib.reload(AstroCell.RGB)
importlib.reload(AstroCell.Image)
importlib.reload(AstroCell.IO)



# Main task
if __name__ == '__main__':

    # Declare whether to operate in parallel (where possible)
    parallel = True

    # State input directory and create output directory inside it
    test_dir = os.path.join(dropbox, 'Work/Scripts/AstroCell/Test/Test_Data/')
    img_dir = 'Histochemial/3100_zeb1/'#'Histochemial/Mammary/Ref_LO_Specific'#'/Flourescant/Liver/APCFLOX1668'#'Histochemial/3100_zeb1/'
    in_dir = os.path.join(test_dir, img_dir)
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

        # Record if operating in parallel
        rgb.RecParallel(parallel)

        # Pass TempDir object to RGB and Image objects
        rgb.TempDir(temp)

        # Clean edges of images
        [ channel.CleanEdges() for channel in rgb.iter ]

        # Create coadd of all three channels
        rgb.MakeCoadd()

        # Create Canny-based edge detection to work out typical cell size
        [ channel.CannyBlobs(sigma=2.0) for channel in rgb.iter_coadd ]

        # Use Canny features in all channels to preliminarily identify regions that hold cells
        rgb.CannyMask()

        # Determine if image is black-background; if not, set it so that it is
        rgb.BlackOnWhite()

        # Use Laplacian-of-Gaussian and Difference-of_Gaussian blob detection to isolate regions occupied by cells
        if parallel:
            logdog_blob_list = joblib.Parallel( n_jobs=mp.cpu_count() )\
                                              ( joblib.delayed( channel.LogDogBlobs )\
                                              ( canny_features=rgb.coadd.canny_features )\
                                              for channel in rgb.iter_coadd )
            [ setattr(rgb.iter_coadd[c],'logdog_mask',logdog_blob_list[c][0]) for c in range(0,len(rgb.iter_coadd)) ]
            [ setattr(rgb.iter_coadd[c],'logdog_features',logdog_blob_list[c][1]) for c in range(0,len(rgb.iter_coadd)) ]
        else:
            [ channel.LogDogBlobs(canny_features=rgb.coadd.canny_features) for channel in rgb.iter_coadd ]

        # Use Canny, LoG, and DoG blobs in all channels to create an improved mask of pixels that represent cells
        rgb.BlobMask()

        # Remove large-scale background structures from image (to create source extraction map)
        rgb.DetFilter()

        # Use canny features to create markers for cells and background, to anchor segmentation
        [ channel.ThreshSegment(rgb.blob_mask) for channel in rgb.iter_coadd ]

        rgb.r.WaterDeblend()

        # Use canny features to create markers for cells and background, to anchor segmentation
        [ channel.WaterDeblend() for channel in rgb.iter_coadd ]



        sdfdfdsvds




# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))
#astropy.io.fits.writeto('/home/chris/coadd.fits', coadd, clobber=True)

# Jubiliate
print('All done!')
