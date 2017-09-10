# Identify location
import socket
location = socket.gethostname()
if location == 'Orthanc':
    dropbox = 'E:\\Users\\Chris\\Dropbox\\'
if location == 'sputnik':
    dropbox = '/home/dhris/Dropbox/'
if location == 'saruman':
    dropbox = '/home/herdata/spx7cjc/Dropbox/'

# Import smorgasbord
import os
import warnings
warnings.simplefilter('ignore', category=Warning)
import matplotlib
matplotlib.use("Pdf")
import matplotlib.pyplot as plt
plt.ioff()
import astropy.logger
astropy.log.setLevel('ERROR')
import AstroCell.Main



# Main process
if __name__ == '__main__':

    # State input directory for test data (various options)
    test_dir = 'Test_Data/'
    #img_dir = 'Histochemial/3100_zeb1/'
    #img_dir = 'Flourescant/Liver/APCFLOX1668/'
    img_dir = 'Flourescant/Mammary/Ref_LO/'
    #img_dir = 'Histochemial/Mammary/Ref_LO/'
    in_dir = os.path.join(test_dir, img_dir)

    # Set output directory for Dills (like pickle jars, these are snapshots to resume AstroCell from a 'saved' point, for testing)
    if location == 'sputnik':
        dill_dir = os.path.join( os.path.expanduser('~'), '/Data/AstroCell/Dills/' )
    else:
        dill_dir = False


    # Launch AstroCell
    AstroCell.Main.Run(in_dir=in_dir, cell_colours=2, substructure_flag=True, parallel=7, mc_factor=1.0, dill_dir=dill_dir, verbose=True)