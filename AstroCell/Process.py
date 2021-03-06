# Import smorgasbord
import os
import pdb
import gc
import copy
import time
import datetime
import warnings
warnings.simplefilter('ignore', category=Warning)
import multiprocessing as mp
import numpy as np
import scipy.stats
import scipy.spatial
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.io.fits
import skimage.feature
import skimage.measure
import webcolors
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



def LabelShuffle(label_map_old, test=False):
    """ Function that takes a labelled segmented map and generates random labels, for more aesthetically pleasing visualisations """

    # Check that no labels are duplicated
    for i in range(1,label_map_old.max()):
        label_map_indv = label_map_old.copy()
        label_map_indv[np.where(label_map_indv!=i)] = 0
        label_map_indv_sublabel = scipy.ndimage.measurements.label(label_map_indv, structure=scipy.ndimage.generate_binary_structure(2,2))
        if label_map_indv_sublabel[1]>1:
            for i in range(2,label_map_indv_sublabel[1]):
                label_map_old[ np.where(label_map_indv_sublabel[0] == i) ] = label_map_old.max() + (i-1)

    # Find each label in map
    label_list = np.unique(label_map_old)
    label_list = label_list[np.where(label_list>0)]
    label_shuffle = np.random.permutation(np.arange(1,label_list.size+1))

    # Loop over labels, picking a new label for each
    label_map_new = np.zeros(label_map_old.shape).astype(int)
    for i in range(1,label_list.shape[0]):
        label_old = label_list[i-1]
        label_new = label_list.shape[0] + label_shuffle[i-1]
        label_map_new[np.where(label_map_old==label_old)] = label_new

    # Shift labels so that they start from 1 again
    if np.where(label_map_new>0)[0].shape[0] > 0:
        label_shift = np.min( label_map_new[ np.where(label_map_new>0) ] ) - 1
        label_map_new[ np.where(label_map_new>0) ] -= label_shift
        label_map_new[ np.where(label_map_new>0) ] += np.nanmax(label_map_new)

    # Return shuffled map
    return label_map_new



def FillHoles(map_in):
    """ A function for robustly filling holes in labelled segments """

    # Loop over labels in segmentation map (skipping null values)
    map_out = map_in.copy().astype(int)
    for i in np.unique(map_in.astype(int)):
        if i<=0:
            continue

        # Fill holes withinin this label only
        map_label = map_in.copy().astype(int)
        map_label[ np.where(map_label!=i) ] = 0
        map_label = scipy.ndimage.morphology.binary_fill_holes(map_label.astype(bool))
        map_out[np.where(map_label==True)] = i

    # Return filled segmentation map
    return map_out



def WaterWrapper(Image, seg_map, iter_total, verbose, img_id):
    """ Wrapper around watershed segmentation function, for ease of parallelisation """

    # Make copy of input Image object to work with
    Image = copy.deepcopy(Image)

    # Decide how many markers to generate, based on number of already-identified features, and the proportion of the map they occupy
    n_thresh_seg = int( np.unique(seg_map).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_canny = int( np.unique(Image.canny_features).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_logdog = int( np.unique(Image.logdog_features).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_markers = int( 3.0 * np.nanmax([n_thresh_seg, n_canny, n_logdog]) )

    # Generate totally random marker coordinates
    markers = np.random.random(size=(n_markers,2))
    markers[:,0] *= Image.map.shape[0]
    markers[:,1] *= Image.map.shape[1]

    # Prune markers located too close to each other
    markers = ProximatePrune(markers, 0.5*np.sqrt(Image.thresh_area/np.pi))

    # Loop over each threshold segment, ensuring that each recieves at least one marker
    thresh_seg_labels = np.unique(seg_map)
    for i in range(1,len(thresh_seg_labels)):
        label = thresh_seg_labels[i]
        label_where = np.where(seg_map==label)
        label_marker = np.random.randint(label_where[0].shape[0])
        markers[i,0] = label_where[0][label_marker]
        markers[i,1] = label_where[1][label_marker]

    # Loop over each Canny Feature segment, ensuring that each recieves at least one marker
    canny_labels = np.unique(Image.canny_features)
    for i in range(0,len(canny_labels)):
        label = canny_labels[i]
        label_where = np.where(Image.canny_features==label)
        label_marker = np.random.randint(label_where[0].shape[0])
        markers[i+len(thresh_seg_labels),0] = label_where[0][label_marker]
        markers[i+len(thresh_seg_labels),1] = label_where[1][label_marker]

    # Convert marker points into a marker array
    marker_map = np.zeros(Image.map.shape, dtype=np.int)
    for i in range(0, markers.shape[0]):
        marker_map[ int(markers[i,0]), int(markers[i,1]) ] = i+1

    # Remove markers that do not lie within segmented objects
    marker_map[np.where(seg_map==0)] = 0

    # Create mask
    mask_map = np.zeros(seg_map.shape).astype(bool)
    mask_map[np.where(seg_map>0)] = True

    """# If input map specified, select it; otherwise, simply use the det map
    if in_map==None:
        in_map = Image.detmap.copy()
    elif (isinstance(in_map, str)) and (in_map in dir(Image)) and (isinstance(getattr(Image, in_map), np.ndarray)):
        in_map = getattr(Image, in_map)
    else:
        raise Exception('Specified in_map is not an attribute of the provided Image object')"""

    # Select input map, and invert it (as watershed requires 'peaks' to become 'valleys')
    in_map = Image.detmap.copy()
    in_map = (-1.0 * in_map) + np.nanmax(in_map)

    # Conduct segmentation
    out_map = skimage.morphology.watershed(in_map, marker_map, connectivity=np.zeros([3,3])+1,
                                         offset=None, mask=None, compactness=0, watershed_line=False)
    out_map[np.where(seg_map==0)] = 0

    # Estimate completion time,
    iter_complete, time_est = ProgressDir(os.path.join(Image.temp_dir,'Prog_Dir'), iter_total, raw=True)

    # Work out when to report estimated completion time, and then do so
    iter_report = np.min([ 10*(mp.cpu_count()-1), int(iter_total/5) ])
    if iter_complete == 1:
        if verbose:
            print('['+img_id+'] Starting Monte-Carlo deblending for '+str(Image.name)+' channel; estimated completion time pending.')
    elif iter_complete == iter_report:
        datetime_now = datetime.datetime.now()
        datetime_est = datetime.datetime.fromtimestamp(float(time_est))
        datetime_delta = datetime_est - datetime_now
        delta_m, delta_s = divmod(datetime_delta.total_seconds(), 60)
        delta_h, delta_m = divmod(delta_m, 60)
        delta_string =  str(int(np.round(delta_h)))+' h, '+str(int(np.round(delta_m)))+' m, '+str(int(np.round(delta_s)))+' s'
        if verbose:
            print('['+img_id+'] Estimated completion time for '+str(Image.name)+' channel Monte-Carlo deblending: '+delta_string+', (at '+str(time.ctime(time_est))+').')#, end='\r')

    # Clean up, and return output segmentation map
    gc.collect
    return out_map



def WalkerWrapper(Image, seg_map, iter_total, verbose, img_id):
    """ Wrapper around random walker segmentation function, for ease of parallelisation """

    # Make copy of input Image object to work with
    Image = copy.deepcopy(Image)

    # Decide how many markers to generate, based on number of already-identified features, and the proportion of the map they occupy
    n_thresh_seg = int( np.unique(seg_map).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_canny = int( np.unique(Image.canny_features).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_logdog = int( np.unique(Image.logdog_features).shape[0] * ( seg_map.size / np.where(seg_map>0)[0].shape[0] ) )
    n_markers = int( 2.0 * np.nanmax([n_thresh_seg, n_canny, n_logdog]) )

    # Generate totally random marker coordinates
    markers = np.random.random(size=(n_markers,2))
    markers[:,0] *= Image.map.shape[0]
    markers[:,1] *= Image.map.shape[1]

    # Prune markers located too close to each other
    markers = ProximatePrune(markers, 0.5*np.sqrt(Image.thresh_area/np.pi))

    # Loop over each threshold segment, ensuring that each recieves at least one marker
    thresh_seg_labels = np.unique(seg_map)
    for i in range(1,len(thresh_seg_labels)):
        label = thresh_seg_labels[i]
        label_where = np.where(seg_map==label)
        label_marker = np.random.randint(label_where[0].shape[0])
        markers[i,0] = label_where[0][label_marker]
        markers[i,1] = label_where[1][label_marker]

    # Loop over each Canny Feature segment, ensuring that each recieves at least one marker
    canny_labels = np.unique(Image.canny_features)
    for i in range(0,len(canny_labels)):
        label = canny_labels[i]
        label_where = np.where(Image.canny_features==label)
        label_marker = np.random.randint(label_where[0].shape[0])
        markers[i+len(thresh_seg_labels),0] = label_where[0][label_marker]
        markers[i+len(thresh_seg_labels),1] = label_where[1][label_marker]

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
    #in_map[np.where(mask_map==False)] = 999

    # Conduct segmentation
    out_map = skimage.segmentation.random_walker(in_map, marker_map, beta=15, mode='cg_mg', tol=0.0025,
                                                      copy=True, multichannel=False, return_full_prob=False, spacing=None)

    # Estimate completion time,
    iter_complete, time_est = ProgressDir(os.path.join(Image.temp_dir,'Prog_Dir'), iter_total, raw=True)

    # Work out when to report estimated completion time, and then do so
    iter_report = np.min([ 10*(mp.cpu_count()-1), int(iter_total/5) ])
    if iter_complete == 1:
        if verbose:
            print('['+img_id+'] Starting Monte-Carlo deblending for '+str(Image.name)+' channel; estimated completion time pending.')
    elif iter_complete == iter_report:
        datetime_now = datetime.datetime.now()
        datetime_est = datetime.datetime.fromtimestamp(float(time_est))
        datetime_delta = datetime_est - datetime_now
        delta_m, delta_s = divmod(datetime_delta.total_seconds(), 60)
        delta_h, delta_m = divmod(delta_m, 60)
        delta_string =  str(int(np.round(delta_h)))+' h, '+str(int(np.round(delta_m)))+' m, '+str(int(np.round(delta_s)))+' s'
        if verbose:
            print('['+img_id+'] Estimated completion time for '+str(Image.name)+' channel Monte-Carlo deblending: '+delta_string+', (at '+str(time.ctime(time_est))+').')#, end='\r')

    # Clean up, and return output segmentation map
    gc.collect
    return out_map



def HysterThresh(in_image, v_low, v_high):
    """ Function (adapted from Emmanuelle Gouillart's blog) to perform Hysteresis thresholding on a probabalistic border map """

    # Create masks of pixels brighter than the high and low hysteresis thresholds
    mask_low = in_image > v_low
    mask_high = in_image > v_high

    # Identift connected pixels
    labels_low = skimage.measure.label(mask_low, background=0) + 1
    count_low = labels_low.max()

    # Check if connected components contain pixels from mask_high
    sums = scipy.ndimage.sum(mask_high, labels_low, np.arange(count_low+1))
    good_label = np.zeros((count_low + 1,), bool)
    good_label[1:] = sums[1:] > 0
    out_mask = good_label[labels_low]

    # Return resulting masl
    return out_mask



def ColourName(requested_colour):
    """ Stack Overflow function by user fraxel to return an English name for an RGB colour (https://tinyurl.com/yav4y9tb) """

    # First see if webcolors has an exact match (unlikely, but hey)
    try:
        closest_name = webcolors.rgb_to_name(tuple(requested_colour))

    # If not exact match, find closest RGB euclidian match
    except ValueError:
        min_colours = {}
        for key, name in webcolors.css3_hex_to_names.items():
            r_c, g_c, b_c = webcolors.hex_to_rgb(key)
            rd = (r_c - requested_colour[0]) ** 2
            gd = (g_c - requested_colour[1]) ** 2
            bd = (b_c - requested_colour[2]) ** 2
            min_colours[(rd + gd + bd)] = name
        closest_name = min_colours[min(min_colours.keys())]

    # Return Ennglish name of best match
    return closest_name

