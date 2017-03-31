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
import AstroCell.Image
import AstroCell.IO
plt.ioff()



# State input directory and create output directory inside it
in_dir = '/home/chris/Data/AstroCell/Histochemial/3100_zeb1/'#'/home/chris/Data/AstroCell/Flourescant/Mammary/Ref_LO'
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
    rgb = AstroCell.Image.RGB(os.path.join(in_dir, in_image))

    # Clean edges of images
    [ channel.CleanEdges() for channel in rgb.iter ]

    # Create coadd of all three channels
    rgb.MakeCoadd()

    # Determine if image is black-background; if not, set it so that it is
    rgb.BlackOnWhite()

    # Remove large-scale background from image (to create source extraction map)
    rgb.DetFilter()

    # Remove *smooth* large-scale background from image (to create flux extraction map)









# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))

# Jubiliate
print('All done!')
