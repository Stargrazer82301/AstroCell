# Import smorgasbord
import pdb
import os
import sys
sys.path.append('/home/chris/Dropbox/Work/Scripts/')
import dill
import imghdr
import numpy as np
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import AstroCell
import AstroCell.RGB
import AstroCell.Image
import AstroCell.IO
plt.ioff()



def Run(in_dir=False, cell_colours=2, substructure_flag=False, parallel=True, mc_factor=1.0, dill_dir=False):
    """ Define function that commences AstroCell cell-counting pipeline """

    # Create output directory
    out_dir = os.path.join(in_dir, 'AstroCell_Output')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        print('[AstroCell] Warning; existing contents of output directory can be overwritten')

    # Initialise parallel status object
    parallel = AstroCell.IO.Parallel(parallel)

    # Identify all image files in input directory
    in_files = os.listdir(in_dir)
    in_files = [in_file for in_file in in_files if not os.path.isdir(os.path.join(in_dir,in_file))]
    in_images = [in_file for in_file in in_files if imghdr.what(os.path.join(in_dir,in_file))!=None]

    # Loop over image files, running them through the pipeline in turn
    for in_image in np.random.permutation(in_images):
        """
        # If testing, load in a pre-processed dill file (to skip uncessary reprocessing)
        rgb = dill.load( open( os.path.join(dill_dir,str('.'.join(in_image.split('.')[:-1]))+'.dj'), 'rb' ) )
        """
        # Initiate AstroCell RGB object
        rgb = AstroCell.RGB.RGB(os.path.join(in_dir,in_image), out_dir)

        # Create temporary directory inside output directory, and store its location
        rgb.TempDir()

        # Record if operating in parallel
        rgb.RecParallel(parallel)

        # Record Monte-Carlo iteration multiplier factor
        rgb.RecMCFactor(mc_factor)

        # Preserve raw, un-modified copies of data for later reference
        [ channel.Raw() for channel in rgb.iter ]

        # Clean edges of images
        [ channel.CleanEdges() for channel in rgb.iter ]

        # Create coadd of all three channels
        rgb.MakeCoadd()

        # Create Canny-based edge detection to work out typical cell size
        [ channel.CannyBlobs(sigma=2.0) for channel in rgb.iter_coadd ]

        # Use Canny features in all channels to preliminarily identify regions that hold cells
        rgb.CannyMask()

        # Determine if image is black-background; if not, set it so that it is
        rgb.BlackOnWhite()

        # Use Laplacian-of-Gaussian and Difference-of-Gaussian blob detection to isolate regions occupied by cells
        rgb.LogDogBlobsWrapper()

        # Use Canny, LoG, and DoG blobs in all channels to create an improved mask of pixels that represent cells
        rgb.BlobMask()

        # Remove large-scale background structures from image (to create source extraction map)
        rgb.DetFilter()

        # Create 'template' cell in each channel, by stacking upon Canny features
        [ channel.CannyCellStack() for channel in rgb.iter_coadd ]

        # Create cross-correlation map in each channel using the 'template' Canny cell
        [ channel.CrossCorr() for channel in rgb.iter_coadd ]

        # Use canny features to create markers for cells and background, to anchor segmentation
        [ channel.ThreshSegment(rgb.blob_mask) for channel in rgb.iter_coadd ]

        # Deblend watershed border maps, to perform segmentations for each band
        [ channel.DeblendSegment() for channel in rgb.iter_coadd ]

        # Combine segments form individual bands to produce final segmentation
        rgb.SegmentCombine()

        # Perform cell 'photometry'
        rgb.CellPhotom()

        # Classify cells
        rgb.CellClassify(cell_colours=cell_colours)

        # Determine actual colours of classified cells
        rgb.CellColours()

        # Create result overview images
        AstroCell.IO.OverviewImages(rgb)

        # Tidy up temporary directory
        rgb.TempDirTidy()

        # Save processed RGB object, for later testing use
        if isinstance(dill_dir, str): rgb.Dill(dill_dir)

    # Report completion
    print('[AstroCell] Processing of all images completed.')





# Main task; generally you want to run AstroCell as a function, but it's useful to initiate a run here way for development and testing
if __name__ == '__main__':

    # Include reloads, to handle any recent changes
    import importlib
    importlib.reload(AstroCell)
    importlib.reload(AstroCell.RGB)
    importlib.reload(AstroCell.Image)
    importlib.reload(AstroCell.IO)

    # State input directories
    test_dir = '/home/chris/Dropbox/Work/Scripts/AstroCell/Test/Test_Data/'
    dill_dir = '/home/chris/Data/AstroCell/Dills/'
    #img_dir = 'Histochemial/3100_zeb1'
    #img_dir = 'Flourescant/Liver/APCFLOX1688_Specific'
    img_dir = 'Flourescant/Mammary/Ref_LO'
    #img_dir = 'Histochemial/Mammary'
    #img_dir = 'Histochemial/Mammary/Ref_LO_Specific'
    in_dir = os.path.join(test_dir, img_dir)

    # Launch AstroCell
    Run(in_dir=False, cell_colours=2, substructure_flag=False, parallel=None, mc_factor=1.0, dill_dir=None)