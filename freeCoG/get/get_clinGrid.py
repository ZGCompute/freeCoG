# get_clinGrid.py
# function to resample to the high-density grid to
# low density (64 contact) clinical recording configuration

import os
import scipy.io
import numpy as np


def get_clinGrid(subj):
    
    # define elecs dir
    elecs_dir = '/Volumes/Thunderbolt_Duo/imaging/subjects/%s/elecs/' %(subj);
    
    # nav to elecs with elecs file
    os.chdir(elecs_dir);
    
    # inds of electrodes recorded in clinical montage
    inds = np.array([[  0,   2,   4,   6,   8,  10,  12,  14,  32,  34,  36,  38,  40,
         42,  44,  46,  64,  66,  68,  70,  72,  74,  76,  78,  96,  98,
        100, 102, 104, 106, 108, 110, 128, 130, 132, 134, 136, 138, 140,
        142, 160, 162, 164, 166, 168, 170, 172, 174, 192, 194, 196, 198,
        200, 202, 204, 206, 224, 226, 228, 230, 232, 234, 236, 238]]);

    # load high density grid
    hd = scipy.io.loadmat('hd_grid.mat');
    hd = hd.get('elecmatrix');

    # resample using indices of electrode coordinates
    clinGrid = hd[inds,:];
    clinGrid = clinGrid.reshape(64,3); #get rid of erroneuos third dimension

    # save resmpled clinical grid
    scipy.io.savemat('clinical_grid.mat', {'elecmatrix':clinGrid});