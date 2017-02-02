# Import smorgasbord
import pdb
import sys
sys.path.append('/home/chris/Dropbox/Work/Scripts/')
import os
import shutil
import imghdr
import inspect
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
import PIL.Image
import AstroCell
import AstroCell.classes
import AstroCell.funcs
plt.ioff()



# State input directory and create output directory inside it
in_dir = '/home/chris/Data/AstroCell/Flourescant/Mammary/Ref_LO/'
out_dir = os.path.join(in_dir, 'AstroCell_Output')
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)
temp = AstroCell.classes.temp_dir(out_dir)
"""
# Create tuple to hold each band object
bands = ( AstroCell.classes.rgb('r'), AstroCell.classes.rgb('g'), AstroCell.classes.rgb('b') )
"""
# Identify and loop over all image files in input directory
in_files = os.listdir(in_dir)
in_files = [in_file for in_file in in_files if not os.path.isdir(os.path.join(in_dir,in_file))]
in_images = [in_file for in_file in in_files if imghdr.what(os.path.join(in_dir,in_file))!=None]
for in_image in in_images:

    # Determine image extension
    in_exten = '.'+in_image.split('.')[-1:][0]

    # Load in bitmap image
    bitmap_image = PIL.Image.open(os.path.join(in_dir, in_image))
    rgb_image = np.array(bitmap_image)

    # Extract each band
    r_image = AstroCell.funcs.edge_clean(rgb_image[:,:,0]).astype(float)
    g_image = AstroCell.funcs.edge_clean(rgb_image[:,:,1]).astype(float)
    b_image = AstroCell.funcs.edge_clean(rgb_image[:,:,2]).astype(float)

    b_edge = skimage.feature.canny(b_image, sigma=1.5).astype(float)

    temp.qw(b_image, in_image.replace(' ','_'))
    temp.qw(b_edge, in_image.replace(' ','_').replace('.png','_canny'))











dfgdfg
# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))

# Jubiliate
print('All done!')