# get_cortSegMeshes.py
# function calls annot2dpv and splitsrf 
# to tessellate and render meshes for freesurfer
# cortical parcellation meshes

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: February 12, 2015

import os
from pialROI_fs2mlab import *

def get_cortSegMeshes(subj, hem):
    
    # location of freesurfer look up table for naming ROIs
    fsLUT_dir= '/Applications/freesurfer/';
    os.chdir(fsLUT_dir);
    
    # open the LUT text file and grab list of label #'s and names
    fname = 'FreeSurferColorLUT.txt';
    file = open(fname, 'r');
    LUT = file.readlines();
    LUT = LUT[429:502]; # grab only crtx ROIs from .txt LUT
    LUT = [row.split("    ") for row in LUT];
    LUT = [row[0:2] for row in LUT]; # only retain label # and name
    LUT.sort();
    LUT.pop(0); # get rid of non-sense string 
    LUT.pop(0); # and get rid of unknown roi
    LUT.pop(35); # rh unknown
    lh_LUT = LUT[0:35];
    rh_LUT = LUT[35:70];
    
    # nav to python scripts dir to import matlab mesh conversion
    python_dir = '/Volumes/Thunderbolt_Duo/Scripts/Python/';
    os.chdir(python_dir);

    
    # nav to subject root dir
    bash_dir = '/Volumes/Thunderbolt_Duo/Scripts/Bash/areal/bin/';
    os.chdir(bash_dir);
    root = '/Volumes/Thunderbolt_Duo/imaging/subjects/' + subj;
    
    # convert annot.aparc to .dpv format
    fs_dir = root + '/fs_data';
    os.system('./annot2dpv %s/label/%s.aparc.annot %s/label/%s.aparc.annot.dpv' %(fs_dir, hem,fs_dir, hem));
    
    # split and render mesh surface for each cortical roi
    os.system('./splitsrf %s/surf/%s.pial.asc %s/label/%s.aparc.annot.dpv %s/surf/%s_pial_roi' %(fs_dir, hem, fs_dir, hem, fs_dir, hem));
    
    # make dir for all pial meshes
    os.chdir((fs_dir + '/surf/'));
    os.mkdir('%s_pial_rois' %(hem));
    os.system('mv %s_pial_roi.* %s_pial_rois/' %(hem, hem));
    os.chdir('%s_pial_rois' %(hem));
    
    # get list of all the .srf roi files created from splitsrf
    srf_list = [i for i in os.listdir(os.getcwd()) if (i.endswith('.srf') and i.startswith(hem))];
    
    if hem == 'lh':
        
       # replace roi number with roi name from LUT
       [os.system('mv %s %s' %(i,(lh_LUT[srf_list.index(i)][1] + '.asc'))) for i in srf_list];
       
    elif hem == 'rh':
        
       # replace roi number with roi name from LUT
       [os.system('mv %s %s' %(i,(rh_LUT[srf_list.index(i)][1] + '.asc'))) for i in srf_list]; # !!!! fix this srf_lsit and rh_LUT do not match
    
    # create corresponding list of ascii file names
    ascii_list = [i for i in os.listdir(os.getcwd()) if (i.endswith('.asc') and i.startswith('ctx-'+ hem))];
    
    # convert ascii to matlab style tri and vert for meshes
    for fname in ascii_list:
       fname = fname;
       pialROI_2mlab(fname);
       
    # mkdirs for tri vert and inds
    os.mkdir('tri')
    os.mkdir('vert')
    os.mkdir('inds')
    os.system('mv *tri.mat tri')
    os.system('mv *vert.mat vert')
    os.system('mv *inds.mat inds')   
    
    
    
    
    

