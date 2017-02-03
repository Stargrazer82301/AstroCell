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
import skimage.restoration
import PIL.Image
import AstroCell
import AstroCell.Image
import AstroCell.IO
plt.ioff()



# State input directory and create output directory inside it
in_dir = '/home/chris/Data/AstroCell/Flourescant/Mammary/Ref_LO'#'/home/chris/Data/AstroCell/Histochemial/3100_zeb1/'
out_dir = os.path.join(in_dir, 'AstroCell_Output')
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)
temp = AstroCell.IO.TempDir(out_dir)



# Identify and loop over all image files in input directory
in_files = os.listdir(in_dir)
in_files = [in_file for in_file in in_files if not os.path.isdir(os.path.join(in_dir,in_file))]
in_images = [in_file for in_file in in_files if imghdr.what(os.path.join(in_dir,in_file))!=None]
for in_image in np.random.permutation(in_images):

    # Determine image extension, and read in image
    exten = '.'+in_image.split('.')[-1:][0]
    bitmap_image = PIL.Image.open(os.path.join(in_dir, in_image))
    rgb_image = np.array(bitmap_image)

    # Load image, and create AstroCell.Image object for each band
    r, g, b = AstroCell.IO.LoadRGB(os.path.join(in_dir,in_image))

    # Determine if image is black-background; if not, set it so that it is

    # Clean edges of images
    [ band.CleanEdges() for band in [r,g,b] ]






    temp.qw(b.map, 'b_image')

    psf = astropy.convolution.Gaussian2DKernel(1.0).array
    b_deconv, _ = skimage.restoration.unsupervised_wiener(b.unclean.map, psf, clip=False)
    temp.qw(b_deconv,'b_deconv')

    b_edge = skimage.feature.canny(b.map, sigma=2.0).astype(float)
    temp.qw(b_edge, 'b_edge')

    b_deconv_edge = skimage.feature.canny(b.map, sigma=2.0).astype(float)
    temp.qw(b_deconv_edge, 'b_deconv_edge')

    b_open1 = scipy.ndimage.grey_opening(b.map, size=(1,1))
    temp.qw(b_open1,'b_open1')

    b_erode1 = scipy.ndimage.grey_erosion(b.unclean.map, size=(1,1))
    temp.qw(b_erode1,'b_erode1')

    close_structure = np.array([[1,1,0],[1,1,1],[1,1,1]])
    b_edge_close = scipy.ndimage.binary_closing(b_edge, structure=close_structure)
    temp.qw(b_edge_close,'b_edge_close')
    b_edge_close_invert = np.invert(b_edge_close)
    temp.qw(b_edge_close_invert,'b_edge_close_invert')



    dsasdsad

    b_edge_dilate = scipy.ndimage.binary_dilation(b_edge)
    b_edge_dilate_fill = scipy.ndimage.binary_fill_holes(b_edge_dilate)
    temp.qw(b_edge_dilate_fill,'b_edge_dilate_fill')


    sdfds















    b_edge_open = scipy.ndimage.binary_opening(b_edge).astype(int)
    b_edge_close = scipy.ndimage.binary_closing(b_edge).astype(int)












# Clean up temporary files
shutil.rmtree(os.path.join(out_dir,'Temp'))

# Jubiliate
print('All done!')
