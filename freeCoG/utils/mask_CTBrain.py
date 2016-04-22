# mask_CTBrain.py
# Script extracts brain mask from high-res T1
# using fsl BET, corregisters the CT to to T1   
# using fsl FLIRT, and applys inverse affine
# to brain mask the CT

import os
from nipype.interfaces.fsl import FLIRT
from nipype.interfaces.fsl import BET
from nipype.interfaces.fsl import ApplyXfm
import nibabel as nib
import cv2
import numpy as np

def mask_CTBrain(CT_dir, T1_dir, frac):
    
    # nav to dir with T1 images
    os.chdir(T1_dir);
    
    # run BET to extract brain from T1
    bet = BET();
    bet.inputs.in_file = T1_dir + '/T1.nii';
    bet.inputs.frac = frac;
    bet.inputs.robust = True;
    bet.inputs.mask = True;
    print "::: Extracting Brain mask from T1 using BET :::";
    bet.run();
    
    # use flrit to correg CT to T1
    flirt = FLIRT();
    flirt.inputs.in_file = CT_dir + '/CT.nii';
    flirt.inputs.reference = T1_dir + '/T1.nii';
    flirt.inputs.cost_func = 'mutualinfo';
    print "::: Estimating corregistration from CT to T1 using FLIRT :::";
    flirt.run();
    
    # get inverse of estimated coreg of CT tp T1
    print "::: Estimating inverse affine for Brain mask of CT :::";
    os.system('convert_xfm -omat inv_CT_flirt.mat -inverse CT_flirt.mat');
    os.system('mv inv_CT_flirt.mat %s' %(CT_dir));
    
    # apply inverse affine coreg to get brain mask in CT space
    applyxfm = ApplyXfm()
    applyxfm.inputs.in_file = T1_dir + '/T1_brain_mask.nii.gz'
    applyxfm.inputs.in_matrix_file =  CT_dir + '/inv_CT_flirt.mat'
    applyxfm.inputs.out_file = CT_dir + '/CT_mask.nii.gz'
    applyxfm.inputs.reference = CT_dir + '/CT.nii';
    applyxfm.inputs.apply_xfm = True
    print "::: Applying inverse affine to Brain mask to get CT mask :::";
    applyxfm.run()
    
    # dilate brain mask to make sure all elecs are in final masked img
    CT_mask = nib.load(CT_dir + '/CT_mask.nii.gz');
    CT_mask_dat = CT_mask.get_data();
    kernel = np.ones((5,5),np.uint8);
    
    print "::: Dilating CT mask :::";
    dilated = cv2.dilate(CT_mask_dat, kernel, iterations = 1);
    hdr = CT_mask.get_header();
    affine = CT_mask.get_affine();
    N = nib.Nifti1Image(dilated,affine,hdr);
    new_fname = CT_dir + '/dilated_CT_mask.nii';
    N.to_filename(new_fname);

    
    # apply mask to CT
    os.chdir(CT_dir);
    print "::: masking CT with dilated brain mask :::";
    os.system("fslmaths 'CT.nii' -mas 'dilated_CT_mask.nii.gz' 'masked_CT.nii'");
    os.system('gunzip masked_CT.nii.gz'); 
    
    
    
    
    
    