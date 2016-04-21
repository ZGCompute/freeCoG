# get_pvc.py
# function takes a nifti volume of functional data,
# and a two grey-matter and white-matter mask
# segmentations extracted from a structural image
# registered w/ the functional volume to calculate 
# a partial volume correction, and normalize the results
# by the mean of WM signal. Binary opening (wm) and closing (gm)
# is performed on the masks to erode noise in the segmentation.

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: June 29, 2015


import os
import nibabel as nib
import numpy as np
from scipy import ndimage

def get_pvc(vol, gm, wm):
    
    
    # load in img data for raw data, 
    img = nib.load(vol);
    img_data = img.get_data();
    img_affine = img.affine;
    img_hdr = img.get_header();
    
    # load gm/wm seg masks
    WM = nib.load(wm);
    WM_data = WM.get_data();
    GM = nib.load(gm);
    GM_data = GM.get_data();
    
    # get rid of nans
    where_nans = np.isnan(img_data);
    img_data[where_nans] = 0;
    where_nans = np.isnan(WM_data);
    WM_data[where_nans] = 0;
    where_nans = np.isnan(GM_data);
    GM_data[where_nans] = 0;
    
    # binary open WM and close GM
    opened_WM = ndimage.binary_opening(WM_data);
    closed_GM = ndimage.binary_closing(GM_data);
    
    # calc partial volume correction
    pvc_data = img_data * ( closed_GM  ) / (  closed_GM + 0.4 * opened_WM + 1e-5 );
    
    # normalize the PVC data by white matter mean signal
    img_wm = img_data* WM_data;
    pvc_norm_data = pvc_data/(img_wm.mean());
    
    # save normalized PVC data to NIFTI volume
    N = nib.Nifti1Image(pvc_norm_data,img_affine,img_hdr);
    fname = 'pvc_wm_norm_' + vol;
    N.to_filename(fname);
    
    