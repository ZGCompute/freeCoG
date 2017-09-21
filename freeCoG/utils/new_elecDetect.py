# new_elecDetect.py
# function uses 3d cv alg to extract electrode center coordinates
# from a post-implant CT, and corregisters the ouputs to a pre-implant
# T1 weighted MRI.

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: April 27, 2015

import os

from CT_electhresh import *
from im_filt import *
from Kmeans_cluster import *
from get_cRAS import *

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
import numpy as np
import cv2
import nibabel as nib
import subprocess
import shlex
from nipype.interfaces import matlab as matlab
from nipype.interfaces import spm
from nipype.interfaces.fsl import FLIRT

def elecDetect(CT_dir, fs_dir,img, thresh_f, thresh_c, frac, num_gridStrips, n_elects):
    
    # run CT_electhresh
    print "::: Thresholding CT at %f stdv :::" %(thresh);
    CT_img = CT_dir + img;
    electhresh(CT_img,thresh_f, thresh_c);
    
    # run im_filt to gaussian smooth thresholded image
    fname = CT_dir + "thresh_" + img;
    print "::: Applying guassian smooth of 2.5mm to thresh-CT :::";
    img_filt(fname);
    
    # compute cv alg to find electrode center coordinates
    mlab = matlab.MatlabCommand();
    mlab.inputs.script = "" %(CT_dir, CT_dir , CT_dir, CT_dir);
    
    print "::: Applying 3D cv alg to smoothed thresh-CT :::";        
    out = mlab.run(); # run the hough and get errors if any ocur

    # save txt version of hough putputs
    coords = scipy.io.loadmat(CT_dir + '/'+ '3dcv_coords.mat').get('coords')
    np.savetxt('dirty_elecsAll.txt', coords, delimiter=' ', fmt='%-7.1f');
    
    # run fsl FLIRT to compute CT-MRI registration
    flirt = FLIRT();
    flirt.inputs.in_file = CT_dir + '/CT.nii';
    flirt.inputs.reference = fs_dir + '/orig.nii';
    flirt.inputs.cost_func = 'mutualinfo';
    print "::: Computing corregistration between CT and MRI using FLIRT :::";
    flirt.run();
            
    # use fsl to calc voxel coords for elec centers in CT registered to orig 
    # in RAS voxel coordinates using affine from FLIRT
    print "::: Applying FLIRT transformation to 3d cv ouput :::";
    orig_brain = fs_dir + '/orig.nii';
    os.system('cat %s/dirty_elecsAll.txt | img2imgcoord -src %s/CT.nii -dest %s \
    -xfm %s/CT_flirt.mat -vox > %s/dirty_elecsAll_RAS.txt' %(CT_dir, CT_dir, orig_brain, fs_dir, CT_dir));
    
    # remove outliers caused by staples in the most superior portion of the img
    print "::: Removing outliers from detected electrode coordinates :::";
    elecs_dirty_RAS = scipy.loadtxt('dirty_elecsAll_RAS.txt', skiprows=1);
    clean_elecs_RAS = [];
    for i in elecs_dirty_RAS:
        if i[1] > (0.3*np.max(elecs_dirty_RAS[:,2])):
            clean_elecs_RAS.append(i);
    clean_elecs_RAS = np.array(clean_elecs_RAS);
    
    # save surface RAS coords for cleaned coords
    # first apply vox2Ras affine for all fs orig volumes
    print "::: Applying freesurfer vox2RAS transformation :::";
    affine = np.array([[  -1.,    0.,    0.,  128.],
           [   0.,    0.,    1., -128.],
           [   0.,   -1.,    0.,  128.],
           [   0.,    0.,    0.,    1.]]);
    
    intercept = np.ones(len(clean_elecs_RAS));
    vox_elecs = np.column_stack((clean_elecs_RAS,intercept));
    vox_elecs = np.transpose(vox_elecs);
    sRAS_coords = np.matrix.dot(affine,vox_elecs);
    sRAS_coords = np.transpose(sRAS_coords);
    sRAS_coords = sRAS_coords[:,0:3];
        
    # apply difference coord (lh_pial.asc  - lh.pial) to shift all elecmatrix coords
    c_RAS = get_cRAS(fs_dir,subj);
    sRAS_coords = coords - c_RAS; # shifted to correct back to lh.pial coords
    scipy.io.savemat('elecs_all_cleaned.mat', {'elecmatrix':sRAS_coords});
    
    # run K means to find grid and strip clusters
    n_clusts = num_gridStrips; # number of cluster to find = number of electrodes
    iters = 1000;
    init_clusts = n_clusts;
    print "::: Performing K-means cluster to find %d grid-strip and noise clusters :::" %(n_clusts);
    [clust_coords, clust_labels] = clust(sRAS_coords, n_clusts, iters, init_clusts);
    
    # save cluster results as individual arrays
    grid_strip_iter = np.linspace(0,num_gridStrips-1, num_gridStrips);
    sRAS_coords_list = [list(b) for b in sRAS_coords];
    clust_labels_list = list(clust_labels);
    for i in grid_strip_iter:
        coords = [];
        for d in sRAS_coords_list:
            if clust_labels_list[sRAS_coords_list.index(d)] == int(i):
               coords.append(d);
        coords = np.array(coords);
        scipy.io.savemat('clust%s.mat' %(i), {'elecmatrix':coords});
        

    #masked_elecs = np.multiply(elecs_img, dilated);
    #clean_elecs_RAS = np.nonzero(masked_elecs);
    #x = np.array(clean_elecs_RAS[0]);
    #y = np.array(clean_elecs_RAS[1]);
    #z = np.array(clean_elecs_RAS[2]);
    #clean_elecs_RAS = np.column_stack((x,y,z));
    
    #save dilated nifti brain mask
    #hdr = brain_mask.get_header();
    #affine = brain_mask.get_affine();
    #N = nib.Nifti1Image(dilated,affine,hdr);
    #new_fname = 'dilated_brain_mask.nii';
    #N.to_filename(new_fname);
    
        # dilate brain mask
    #kernel = np.ones((12,12),np.uint8);
    #brain_mask = nib.load('brain_mask.nii');
    #brain_mask_dat = brain_mask.get_data();
    #dilated = cv2.dilate(brain_mask_dat, kernel, iterations =1);
    
    # Extract and apply CT brain mask to hough voxel output: need to edit this
    #orig_brain = fs_dir + 'mri/orig.nii';
    #os.system("thresh_it=`fslstats '%s' -P 15`" %(orig_brain));
    #os.system("'fslmaths 'brain.nii' -thr $thresh_it -bin brain_mask.nii");
