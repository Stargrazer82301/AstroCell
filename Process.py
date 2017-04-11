# Import smorgasbord
import pdb
import time
import numpy as np
import scipy.stats
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.io.fits
import skimage.feature
from ChrisFuncs import SigmaClip
plt.ioff()





def CannyCells(in_image, sigma=2.0):
    """ A function that uses Canny edge detection to generate a initial guess of what regions of an image are occupied by cells """

    # Evaluate noise in image, via iterative sigma-clipping
    clip = SigmaClip(in_image, median=True, sigma_thresh=3.0)

    # Use noise level to generate atificial noise to add to image
    noise = np.random.normal(0.0, 0.5*clip[0], size=(in_image.shape[0],in_image.shape[1],60))

    # Add each generation of random noise to image in turn, performing Canny edge filtering on each
    canny_stack = np.zeros(noise.shape)
    for i in range(0,noise.shape[2]):
        canny_stack[:,:,i] = skimage.feature.canny(in_image+noise[:,:,i], sigma=sigma)

    # Co-add each individual Canny iteration; Canny edges found in lots of iterations are assumed to be real
    canny = np.sum(canny_stack, axis=2)
    canny[np.where(canny<20)] = 0
    canny[np.where(canny>=20)] = 1
    #astropy.io.fits.writeto('/home/chris/canny.fits', canny.astype(float), clobber=True)

    # Run map through one iteration of binary closing, to close up gaps in the perimeters of canny edges
    canny_close = scipy.ndimage.binary_closing(canny, iterations=1, structure=scipy.ndimage.generate_binary_structure(2,1))

    # Label closed Canny image
    canny_cell_labels = skimage.measure.label(np.invert(canny_close), connectivity=1)

    # Record number of pixels in each labelled feature
    canny_cell_areas = np.unique(canny_cell_labels, return_counts=True)[1]
    canny_cell_areas = canny_cell_areas[np.where(canny_cell_areas>0)]

    # Clip label areas, to identify threshold above which cell regions are large enough to likely be spurious
    labels_clip = SigmaClip(canny_cell_areas, median=True, sigma_thresh=4.0)
    labels_area_thresh = np.max(labels_clip[2])
    labels_exclude = np.arange(0,canny_cell_areas.size)[ np.where( (canny_cell_areas>labels_area_thresh) | (canny_cell_areas<=5) ) ]

    # Remove spurious labels (flattening to improve speed, then reshaping after processing)
    canny_cells = canny_cell_labels.copy().flatten()
    canny_cells[np.in1d(canny_cells,labels_exclude)] = 0
    canny_cells = np.reshape(canny_cells, canny_cell_labels.shape)

    # Return final image
    return canny_cells



def SegmentCombine():
    """ A function that combines segmentation maps from multiple channels to create a master segmentation map """


