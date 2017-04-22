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
    n_thresh_seg = int( np.unique(seg_map).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_canny = int( np.unique(Image.canny_features).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_logdog = int( np.unique(Image.logdog_features).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_markers = 2.0 * np.nanmax([n_thresh_seg, n_canny, n_logdog])

    # Generate totally random marker coordinates
    markers = np.random.random(size=(n_markers,2))
    markers[:,0] *= Image.map.shape[0]
    markers[:,1] *= Image.map.shape[1]

    # Prune markers located too close to each other
    markers = ProximatePrune(markers, 0.5*np.sqrt(Image.thresh_area/np.pi))

    # Loop over each segment, ensuring that each recieves at least one marker
    thresh_seg_labels = np.unique(seg_map)
    for i in range(1,len(thresh_seg_labels)):
        label = thresh_seg_labels[i]
        label_where = np.where(seg_map==label)
        label_marker = np.random.randint(label_where[0].shape[0])
        markers[i,0] = label_where[0][label_marker]
        markers[i,1] = label_where[1][label_marker]

    # Convert marker points into a marker array
    marker_map = np.zeros(Image.map.shape, dtype=np.int)
    for i in range(0, markers.shape[0]):
        marker_map[ int(markers[i,0]), int(markers[i,1]) ] = i+1

    # Remove markers that do not lie within segmented objects
    marker_map[np.where(seg_map==0)] = 0

    # Create mask and invert map
    mask_map = np.zeros(seg_map.shape).astype(bool)
    mask_map[np.where(seg_map>0)] = True
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


