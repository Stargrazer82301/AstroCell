# Import smorgasbord
import sys
import os
import warnings
warnings.simplefilter('ignore', category=Warning)
import numpy as np
import matplotlib.pyplot as plt
import astropy.logger
astropy.log.setLevel('ERROR')
plt.ioff()
import AstroCell.Main



# Main process
if __name__ == '__main__':

    # State input directory for test data (various options)
    test_dir = 'Test_Data/'
    #img_dir = 'Histochemial/3100_zeb1'
    #img_dir = 'Flourescant/Liver/APCFLOX1668/'
    #img_dir = 'Flourescant/Mammary/Ref_LO/Specific'
    img_dir = 'Histochemial/Mammary/Ref_LO'
    in_dir = os.path.join(test_dir, img_dir)

    # Set output directory for Dills (like pickle jars, these are snapshots to resume AstroCell from a 'saved' point, for testing)
    dill_dir = False
    dill_dir = os.path.join( os.path.expanduser('~'), '/Data/AstroCell/Dills/' )

    # Launch AstroCell
    AstroCell.Main.Run(in_dir=in_dir, cell_colours=2, substructure_flag=False, parallel=4, mc_factor=1.0, dill_dir=dill_dir, verbose=True)