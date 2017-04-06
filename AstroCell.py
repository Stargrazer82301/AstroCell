# Import smorgasbord
import sys
import os
import shutil
import PIL.Image
import imghdr
import astropy.convolution
import astropy.stats
import astropy.visualization
import astropy.visualization.mpl_normalize
import astropy.io.fits
import numpy as np
import matplotlib.pylab as plt
import photutils
import AstroCell.Image
import AstroCell.IO
plt.ioff()





def Count(in_dir):
    """ Define function that commences AstroCell cell-counting pipeline """



    # Create output directory inside target directory
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
    for in_image in in_images:

        # Read in raw image, constructing an AstroCell RGB object
        rgb = AstroCell.RGB.RGB(os.path.join(in_dir, in_image))

        # Clean edges of images
        [ channel.CleanEdges() for channel in rgb.iter ]

        # Create coadd of all three channels
        rgb.MakeCoadd()

        # Determine if image is black-background; if not, set it so that it is
        rgb.BlackOnWhite()

        # Remove large-scale background structures from image (to create source extraction map)
        rgb.DetFilter()

        # Construct basic matched filter in each channel, using Canny features
        [ channel.CannyCellStack() for channel in rgb.iter_coadd ]

        # Use canny features to create markers for cells and background, to anchor segmentation
        [ channel.ThreshSegment(rgb.canny_mask) for channel in rgb.iter_coadd ]















