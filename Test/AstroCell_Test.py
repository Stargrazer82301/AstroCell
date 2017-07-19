# Import smorgasbord
import pdb
import sys
sys.path.append('/home/chris/Dropbox/Work/Scripts/')
import os
import shutil
import inspect
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pylab as plt
import astropy.logger
astropy.log.setLevel('ERROR')
import astropy.convolution
import astropy.stats
import astropy.visualization
import astropy.visualization.mpl_normalize
import astropy.io.fits
import AstroCell
import AstroCell.Main
import AstroCell.RGB
import AstroCell.Image
import AstroCell.IO
plt.ioff()

# Identify location and set Dropbox path accordingly
import socket
location = socket.gethostname()
if location == 'Monolith':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'sputnik':
    dropbox = '/home/chris/Dropbox/'
if location == 'saruman':
    dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Include reloads, to handle any recent changes
import importlib
importlib.reload(AstroCell)
importlib.reload(AstroCell.RGB)
importlib.reload(AstroCell.Image)
importlib.reload(AstroCell.IO)
importlib.reload(AstroCell.Process)



# Main process
if __name__ == '__main__':

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
    AstroCell.Main.Run(in_dir=in_dir, cell_colours=2, substructure_flag=False, parallel=True, mc_factor=1.0, dill_dir=dill_dir)