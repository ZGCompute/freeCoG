# get_cRAS.py

# function gets freesurfer c_RAS vector
# to save for shifting meshes into orig space correctly


import os
import nibabel as nib
import subprocess
import shlex
import scipy.io
import numpy as np

def get_cRAS(fs_dir,subj):
    
    # define subject dir
    subj_fs_dir = fs_dir + '/';
    
    # use mris_info to get c_RAS from lh.pial freesurfer ouput
    os.system('mris_info %s/surf/lh.pial >& %s/surf/lh_pial_info.txt' %(subj_fs_dir, subj_fs_dir));
    command = 'grep "c_(ras)" %s/surf/lh_pial_info.txt' %(subj_fs_dir);
    args = shlex.split(command);
    p = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0];
    p = p.strip('\n'); p = p.strip('c_(ras): '); p = p.strip('(');p = p.strip(')');
    p = p.split(','); p = [float(i) for i in p]; p = np.array(p);
    
    c_RAS = p; # this is the coordinate shift from lh.pial to lh.pial.asc coords
    scipy.io.savemat('c_RAS.mat', {'c_RAS':c_RAS});
    
    return c_RAS;