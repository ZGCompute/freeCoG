# make_fs_CortROIs.py
# function takes a patient's freesurfer aparc+aseg.mgz 
# cortical and subcortical segmentation labeled volume
# and produces individual volumetric ROIs for each 
# region to be used for region-region connectomics

# Author: Zachary Greenberg
# UCSF Department of Neurosurgery
# Date Created: 11/26/2014

import os
import nipype
import subprocess
import shlex

def make_CortROIs(subj):
    
    # location of freesurfer look up table for naming ROIs
    fsLUT_dir= '/Applications/freesurfer/';
    os.chdir(fsLUT_dir);
    
    # open the LUT text file and grab list of label #'s and names
    fname = 'FreeSurferColorLUT.txt';
    file = open(fname, 'r');
    LUT = file.readlines();
    cort_LUT = LUT[429:502]; # grab only crtx ROIs from .txt LUT
    cort_LUT = [row.split("    ") for row in cort_LUT];
    cort_LUT = [row[0:2] for row in cort_LUT]; # only retain label # and name
    cort_LUT.sort();
    cort_LUT.pop(0); # get rid of non-sense string
    
    # get subcort seg lists too
    lh_subcort = [['8', ' Left-Cerebellum-Cortex'],
    ['9', ' Left-Thalamus'],
    ['10', 'Left-Thalamus-Proper'],
    ['11', 'Left-Caudate'],
    ['12', 'Left-Putamen'],
    ['13', 'Left-Pallidum'],
    ['16', 'Brain-Stem'],
    ['17', 'Left-Hippocampus'],
    ['18', 'Left-Amygdala'],
    ['19', 'Left-Insula'],
    ['20', 'Left-Operculum'],
    ['26', 'Left-Accumbens-area'],
    ['27', 'Left-Substancia-Nigra'],
    ['28', 'Left-VentralDC']];
    
    # r hem subcort
    rh_subcort = [['47', 'Right-Cerebellum-Cortex'],
    ['48', 'Right-Thalamus'],
    ['49', 'Right-Thalamus-Proper'],
    ['50', 'Right-Caudate'],
    ['51', 'Right-Putamen'],
    ['52', 'Right-Pallidum'],
    ['53', 'Right-Hippocampus'],
    ['54', 'Right-Amygdala'],
    ['55', 'Right-Insula'],
    ['56', 'Right-Operculum'],
    ['58', 'Right-Accumbens-area'],
    ['59', 'Right-Substancia-Nigra'],
    ['60', 'Right-VentralDC']];
    
    # once list with all roi #s and names
    roi_list = cort_LUT + lh_subcort + rh_subcort;

    # nav to patient's fs_data/mri dir 
    imgs_dir = '/Volumes/Thunderbolt_Duo/imaging/subjects/';
    patient_dir = imgs_dir + subj + '/fs_data/mri';
    os.chdir(patient_dir);
    
    # convert aparc+aseg to NIFTI
    os.system('mri_convert aparc+aseg.mgz aparc_aseg.nii');
    
    # loop through LUT to make roi for each label
    for row in roi_list:
        label_num = int(row[0]);
        label_name = row[1];
        print "::: Creating ROI for %s :::" %(label_name);
        command1 = 'fslmaths aparc_aseg.nii -uthr %d -thr %d %s.nii' %(label_num, label_num, label_name);# create label ROI
        args1 = shlex.split(command1);
        p1 = subprocess.Popen(args1);
        command2 = 'fslmaths %s.nii.gz -div %d %s.nii.gz' %(label_name, label_num, label_name); # binarize mask to 1s and 0s
        args2 = shlex.split(command2);
        p2 = subprocess.Popen(args1);
        
    # make dir to move all labels into
    os.mkdir('crtx_subcort_ROIs');
    os.system('mv *.nii.gz crtx_subcort_ROIs/');
    
    print "::: Finished :::";
    
    