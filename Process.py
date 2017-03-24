# Import smorgasbord
import pdb
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

    # Pass map through canny filter
    canny = skimage.feature.canny(in_image, sigma=1.0)

    # Run map through one iteration of binary closing, to close up gaps in the perimeters of canny edges
    canny_close = scipy.ndimage.binary_closing(canny)

    # Label closed Canny image
    canny_close_labels = skimage.measure.label(np.invert(canny_close), connectivity=1)

    # Loop over labels, recording number of pixels in each
    labels_areas = np.zeros(np.max(canny_close_labels))
    for i in range(0, np.max(canny_close_labels)):
        labels_areas[i] = np.where(canny_close_labels==i)[0].size

    # Clip label areas, to identify threshold above which cell regions are large enough to likely be spurious
    labels_clip = SigmaClip(labels_areas, median=True, sigma_thresh=4.0)
    labels_area_thresh = np.max(labels_clip[2])
    labels_exclude = np.arange(0,labels_areas.size)[ np.where( (labels_areas>labels_area_thresh) | (labels_areas<=5) ) ]

    # Remove spurious labels
    canny_cells = canny_close_labels.copy()
    for label in labels_exclude:
        canny_cells[ np.where( canny_close_labels == label ) ] = 0

    # Return final image
    return canny_cells


