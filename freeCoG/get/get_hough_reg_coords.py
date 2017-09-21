# get_reg_coords.py

import os
import numpy as np
import nibabel as nib
import scipy.io

def reg_coords(coords_img,source, target):
    
    # load in detected electrode coordinates
    # from source CT img
    coords_img = nib.load(coords_img);
    coords_data = coords_img.get_data();
    
    # grab 1's from img representing detected 
    # electrode centers from cv
    coords = np.nonzero(coords_data);
    x = np.array(coords[0]);
    y = np.array(coords[1]);
    z = np.array(coords[2]);
    vox_coords = np.column_stack((x,y,z));
    
    # save electrode centers as voxel coords in native CT space
    scipy.io.savemat('cv_vox_coords.mat', {'elecmatrix':vox_coords});
    
    # load in CT2MM affine matrix
    CT_img = nib.load(source);
    CT_data = CT_img.get_data();
    CT2mm = CT_img.affine;
    
    # load in target MRI affine
    MRI_img = nib.load(target);
    MRI_data = MRI_img.get_data();
    MRI_hdr = MRI_img.get_header();
    MRI2mm = MRI_img.affine;
    CT2MRI = np.linalg.inv(MRI2mm).dot(CT2mm);
    
    # apply CT2MRI affine registration
    reg_coords = nib.affines.apply_affine(CT2MRI,coords);
    scipy.io.savemat('cv_reg_coords.mat', {'elecmatrix':reg_coords});
    
    # save hough reg voxel coords in a NIFTI img for QC
    elecs_img = np.zeros(MRI_data.shape);
    for i in reg_coords:
       elecs_img[i[0], i[1], i[2]] = 1;
    N = nib.Nifti1Image(elecs_img,MRI2mm,MRI_hdr);
    N.to_filename('cv_reg_coords.nii.gz');
    
    # vox2Ras affine for all fs orig volumes
    vox2RAS = np.array([[  -1.,    0.,    0.,  128.],
          [   0.,    0.,    1., -128.],
          [   0.,   -1.,    0.,  128.],
          [   0.,    0.,    0.,    1.]]);
    
    # apply vox2RAS to get surface mesh (sRAS) coordinates      
    sRAS_coords = nib.affines.apply_affine(vox2RAS,reg_coords);
    scipy.io.savemat('cv_sRAS_coords.mat', {'elecmatrix':sRAS_coords});
