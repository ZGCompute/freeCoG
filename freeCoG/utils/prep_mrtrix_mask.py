# prep_mrtrix_mask.py

# script takes dtiInit output DWI and produces masked DWI and binary mask

from dipy.segment.mask import median_otsu
import numpy as np 
import nibabel 
import os

def get_DWImask(DWI_dir):
    
   #nav to DWI dir
   os.chdir(DWI_dir);
    
   # read in DWI
   img = nibabel.load('dwi_aligned_trilin.nii.gz'); 
   data = img.get_data(); 

   # create mask with dipy
   dwi_mask, mask = median_otsu(data, 2, 1);
   mask_img = nibabel.Nifti1Image(mask.astype(np.float32), img.get_affine());
   dwi_img = nibabel.Nifti1Image(dwi_mask.astype(np.float32), img.get_affine());

   # save mask and dwi masked
   nibabel.save(mask_img, 'binary_mask.nii.gz');
   nibabel.save(dwi_img, 'dwi_mask.nii.gz');
