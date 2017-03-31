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





def CannyCells(in_image, sigma=1.0):
    """ A function that uses Canny edge detection to generate a initial guess of what regions of an image are occupied by cells """

    # Pass image through a canny filter
    canny = skimage.feature.canny(in_image, sigma=sigma)

    # Run map through one iteration of binary closing, to close up gaps in the perimeters of canny edges
    canny_close = scipy.ndimage.binary_closing(canny)

    # Label closed Canny image
    canny_close_labels = skimage.measure.label(np.invert(canny_close), connectivity=1)

    # Loop over labels, recording number of pixels in each
    canny_close_labels_flat = canny_close_labels.flatten()
    labels_areas = np.zeros(np.max(canny_close_labels_flat))
    for i in range(0, np.max(canny_close_labels_flat)):
        labels_areas[i] = np.where(canny_close_labels_flat==i)[0].size

    # Clip label areas, to identify threshold above which cell regions are large enough to likely be spurious
    labels_clip = SigmaClip(labels_areas, median=True, sigma_thresh=4.0)
    labels_area_thresh = np.max(labels_clip[2])
    labels_exclude = np.arange(0,labels_areas.size)[ np.where( (labels_areas>labels_area_thresh) | (labels_areas<=5) ) ]

    # Remove spurious labels (flattening to improve speed, then reshaping after processing)
    canny_cells = canny_close_labels.copy().flatten()
    canny_cells[np.in1d(canny_cells,labels_exclude)] = 0
    canny_cells = np.reshape(canny_cells, canny_close_labels.shape)

    # Return final image
    return canny_cells


