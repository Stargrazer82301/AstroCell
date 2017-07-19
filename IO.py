# Import smorgasbord
import pdb
import sys
import os
import shutil
import multiprocessing as mp
import numpy as np
import scipy.stats
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.io.fits
import skimage.segmentation
import PIL.Image
import AstroCell.Image
plt.ioff()





def LoadRGB(in_path):
    """ Function that reads in an image file and returns three AstrCell.Image objects (r, g, b) """

    # Read in image, and conver to array
    bitmap_image = PIL.Image.open(in_path)
    rgb_image = np.array(bitmap_image)

    # Create an Image object from each channel in turn
    r = AstroCell.Image.Image(rgb_image[:,:,0])
    g = AstroCell.Image.Image(rgb_image[:,:,1])
    b = AstroCell.Image.Image(rgb_image[:,:,2])

    # Return Image objects
    return r, g, b



class Parallel():
    """ Class that stores whether, and what kind, of parallel operation is desired """



    def __init__(self, parallel):
        """ Initialise the Parallel object """

        # Parse and record parallel request information
        if parallel == False:
            self.parallel = False
            self.threads = 1
            self.subthreads = 1
        elif isinstance(parallel,float) or isinstance(parallel,float):
            self.threads = 1
            self.parallel = True
            self.subthreads = int(parallel)
        elif parallel == None:
            self.threads = 1
            self.parallel = True
            self.subthreads = mp.cpu_count()-1



def OverviewImages(RGB):
    """ Method that produces overview image figure of the results """

    # Use size of image, and number of panes desired, to determine size of figure
    map_aspect = float(RGB.coadd.map.shape[1]) / float(RGB.coadd.map.shape[0])
    fig_x_panes = 1
    fig_y_panes = 2
    fig_aspect = ( fig_x_panes * map_aspect ) / fig_y_panes
    fig_x_dim = 10.0 * fig_aspect
    fig_y_dim = 10.0

    # Create figure and axes
    fig, axes = plt.subplots(fig_y_panes, fig_x_panes, figsize=(fig_x_dim, fig_y_dim))

    # Remove ticks and frames from axes
    [ ax.set_xticklabels([]) for ax in axes ]
    [ ax.set_yticklabels([]) for ax in axes ]
    [ ax.axhline(linewidth=5, color='black') for ax in axes ]
    [ ax.axvline(linewidth=5, color='black') for ax in axes ]
    [ ax.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off') for ax in axes ]

    # Generate dilated verison of seg map, for display
    image_label_dilate = np.zeros(RGB.segmap.shape) - 1
    for i in np.unique(RGB.labels_list):
        for j in range(0, len(RGB.table)):
            if RGB.table['label'][j] == i:
                image_label_indv = np.zeros(RGB.segmap.shape)
                image_label_indv[ np.where(RGB.segmap==RGB.table['id'][j]) ] = 1
                image_label_indv = scipy.ndimage.morphology.binary_erosion(image_label_indv)
                image_label_dilate[ np.where(image_label_indv==1) ] = i
    image_label_dilate[np.where(RGB.segmap==0)] = -1

    # Plot label borders onto image, for each label
    image_label_border = RGB.cube.copy() / 255
    for i in np.unique(RGB.labels_list):
        image_label_temp = image_label_dilate.copy()
        image_label_temp[np.where(image_label_temp!=i)] = -1
        image_label_border = skimage.segmentation.mark_boundaries(image_label_border, image_label_temp,
                                                                  color=RGB.labels_rgb_bright[i,:], mode='thick', background_label=-1)

    # Shade in greyscale image regions according to label
    image_label_colour = skimage.color.label2rgb(RGB.labelmap, image=image_label_border,
                                                 colors=RGB.labels_rgb_bright.tolist(), bg_label=-1, bg_color=[1,1,1], image_alpha=0.999)

    # Plot image panels
    axes[0].imshow(image_label_border, origin='upper')
    axes[1].imshow(image_label_colour, origin='upper')

    # Set figure to have tight layout, remove axis frames, and save to file
    fig.tight_layout()
    fig.savefig( os.path.join( '/home/chris/', '.'.join(RGB.in_file.split('.')[:-1])+'_output.png' ), dpi=200 ) #RGB.out_dir
    #astropy.io.fits.writeto('/home/chris/bob.fits', image_label_dilate.astype(float), clobber=True)



class Bunch:
    """ Convenience class for gathering related data in an object """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)





