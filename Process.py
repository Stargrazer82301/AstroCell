# Import smorgasbord
import os
import pdb
import gc
import copy
import numpy as np
import scipy.stats
import scipy.spatial
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.io.fits
import skimage.feature
from ChrisFuncs import SigmaClip, ProgressDir
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

    # Label closed Canny image features
    canny_cell_labels = skimage.measure.label(np.invert(canny_close), connectivity=1)

    # Record number of pixels in each labelled feature
    canny_cell_areas = np.unique(canny_cell_labels, return_counts=True)[1]
    canny_cell_areas = canny_cell_areas[np.where(canny_cell_areas>0)]

    # Clip label areas, to identify threshold above which cell regions are large enough to likely be spurious
    labels_clip = SigmaClip(canny_cell_areas, median=True, sigma_thresh=5.0)
    labels_area_thresh = np.max(labels_clip[2])
    labels_exclude = np.arange(0,canny_cell_areas.size)[ np.where( (canny_cell_areas>labels_area_thresh) | (canny_cell_areas<=5) ) ]

    # Remove spurious labels (flattening to improve speed, then reshaping after processing)
    canny_cells = canny_cell_labels.copy().flatten()
    canny_cells[np.in1d(canny_cells,labels_exclude)] = 0
    canny_cells = np.reshape(canny_cells, canny_cell_labels.shape)

    # Fill cell regions, dilate them by 1 pixel (to account for the width of the Canny border), then relabel
    canny_fill = scipy.ndimage.binary_fill_holes(canny_cells, structure=scipy.ndimage.generate_binary_structure(2,1))
    canny_dilate = scipy.ndimage.binary_dilation(canny_fill, iterations=1, structure=scipy.ndimage.generate_binary_structure(2,2))
    canny_cells = skimage.measure.label(canny_dilate, connectivity=1)
    canny_cells = LabelShuffle(canny_cells)
    #astropy.io.fits.writeto('/home/chris/canny_cells.fits', canny_cells, clobber=True)

    # Catch and fix when most of map is a 'feature'
    mode = scipy.stats.mode(canny_cells.flatten())[0][0]
    if mode != 0:
        mode_frac = np.where(canny_cells==mode)[0].shape[0] / canny_cells.size
        if mode_frac > 0.5:
            canny_cells = np.zeros(canny_cells.shape).astype(int)

    # Return final image
    return canny_cells



def ProximatePrune(points, distance):
    """ Function that takes an array containing 2D coordinats, and removes any that lie within a given distance of other points using a KD tree """

    # Construct KD tree, and find pairs within exclusion distance
    tree = scipy.spatial.KDTree(points)
    pairs = tree.query_pairs(distance)

    # Loop over pairs, pruning the first member encountered
    prune = []
    for pair in pairs:
        if pair[0] in prune:
            continue
        if pair[1] in prune:
            continue
        else:
            prune.append(pair[1])

    # Return un-pruned points
    keep = np.setdiff1d( np.arange(points.shape[0]), prune )
    points = np.array([ points[:,0][keep], points[:,1][keep] ]).transpose()
    return points



def LabelShuffle(label_map_old):
    """ Function that takes a labelled segmented map and generates random labels, for more aesthetically pleasing visualisations """

    # Find each label in map
    label_list = np.unique(label_map_old)
    label_shuffle = np.random.permutation(label_list[np.where(label_list>0)])
    label_map_new = label_map_old.copy()

    # Loop over labels, picking a new label for each
    for label_old in label_list:
        if label_old==0:
            continue
        label_new = label_shuffle.shape[0] + np.random.choice(label_shuffle, replace=False)
        label_map_new[np.where(label_map_old==label_old)] = label_new

    # Return shuffled map
    return label_map_new



def WaterWrapper(Image, seg_map, iter_total):
    """ Wrapper around watershed segmentation function, for ease of parallelisation """

    # Make copy of input Image object to work with
    Image = copy.deepcopy(Image)

    # Decide how many markers to generate, based on number of already-identified features, and the proportion of the map they occupy
    n_markers = 1000#n_markers = int( 1.0 * np.unique(seg_map).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] )

    # Generate marker coordinates
    markers = np.random.random(size=(n_markers,2))
    markers[:,0] *= Image.map.shape[0]
    markers[:,1] *= Image.map.shape[1]

    # Prune markers located too close to each other
    markers = ProximatePrune(markers, 0.5*np.sqrt(Image.thresh_area/np.pi))

    # Convert marker points into a marker array
    marker_map = np.zeros(Image.map.shape, dtype=np.int)
    for i in range(0, markers.shape[0]):
        marker_map[ int(markers[i,0]), int(markers[i,1]) ] = i+1

    """# Remove markers that do not lie within segmented objects
    marker_map[np.where(seg_map==0)] = 0"""

    # Create mask and invert map and then conduct segmentation
    mask_map = np.zeros(seg_map.shape).astype(bool)
    mask_map[np.where(seg_map>0)] = True


    mask_map = None


    in_map = Image.detmap.copy()
    in_map = (-1.0 * in_map) + np.nanmax(in_map)

    # Conduct segmentation
    out_map = skimage.morphology.watershed(in_map, marker_map, connectivity=1,
                                         offset=None, mask=mask_map, compactness=0, watershed_line=False)

    # Estimate completion time
    iter_complete, time_est = ProgressDir(os.path.join(Image.temp.dir,'Prog_Dir'), iter_total)
    print('Monte-Carlo deblending iteration '+str(iter_complete)+' of '+str(iter_total)+'; estimated completion time '+str(time_est)+'.')

    # Clean up, and return output segmentation map
    gc.collect
    return out_map



def SegmentCombine():
    """ A function that combines segmentation maps from multiple channels to create a master segmentation map """


