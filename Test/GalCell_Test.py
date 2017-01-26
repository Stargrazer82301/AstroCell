# Import smorgasbord
import pdb
import sys
import os
import shutil
import imghdr
import inspect
import numpy as np
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
plt.ioff()



class temp_dir():
    """ Convenience class for working with temporary directory"""

    def __init__(self, out_dir):
        """ Initialise the object """
        self.dir = os.path.join(out_dir,'Temp')
        if not os.path.exists(self.dir):
            os.mkdir(self.dir)

    def qw(self, data, name):
        """ Method for quickly writing a numpy array to a FITS file in the temporary directory """
        if name[-4:]!='.fits':
            if name[-1:]=='.':
                name = name[:-1] + '.fits'
            else:
                name += '.fits'
        astropy.io.fits.writeto(os.path.join(self.dir,name), data, clobber=True)



# State input directory and create output directory inside it
in_dir = '/home/chris/Data/GalCell/Flourescant/Mammary/Ref_LO/'
out_dir = os.path.join(in_dir, 'GalCell_Output')
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)
temp = temp_dir(out_dir)

# Identify and loop over all image files in input directory
in_files = os.listdir(in_dir)
in_files = [in_file for in_file in in_files if not os.path.isdir(os.path.join(in_dir,in_file))]
in_images = [in_file for in_file in in_files if imghdr.what(os.path.join(in_dir,in_file))!=None]
for in_image in in_images:

    # Determine image extension
    in_exten = '.'+in_image.split('.')[-1:][0]

    # Load in bitmap image an convert to arrays
    bitmap_image = PIL.Image.open(os.path.join(in_dir, in_image))
    rgb_image = np.array(bitmap_image)
    r_image = rgb_image[:,:,0].astype(float)
    g_image = rgb_image[:,:,1].astype(float)
    b_image = rgb_image[:,:,2].astype(float)

    b_edge = skimage.feature.canny(b_image, sigma=2.0).astype(float)
    sdfdsfds







# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))

# Jubiliate
print('All done!')
