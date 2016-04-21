# get_mrtrix_metrics.py
# Script extracts tensor models, FA, and major eigenvectors from HARDI diffusion 
# data pre-processed with vistasoft-AFQ dtiInit pipeline (disotrtion-eddy current 
# correction), and then extracts a response function from the DWI, and uses it to 
# perfor constrained spherical deconvolution and probobailistic fiber tracking
# on the resultant FOD fields

# Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Created: MArch 13, 2015

import os
os.chdir('/Volumes/Thunderbolt_Duo/Scripts/Python/');
from prep_mrtrix_mask import *

def get_mrtrix_metrics(DWI_dir):
    
    # nav to DWI dir
    os.chdir(DWI_dir);
    
    # use dipy to pull out white matter mask
    get_DWImask(DWI_dir);
    
    # convert wm masked dwi to mrtrix .mif format
    os.system('mrconvert dwi_mask.nii.gz dwi.mif');
 
    # convert binarized mask to .mif too
    os.system('mrconvert binary_mask.nii.gz mask.mif');
    
    # compute tensor model from dwi HARDI data
    os.system('dwi2tensor dwi.mif dt.mif -mask mask.mif -grad grads');
    
    # compute FA map from tensor model
    os.system('tensor2metric dt.mif -fa fa.mif -mask mask.mif');
    
    # compute eigenvector rgb map from tensor
    os.system('tensor2metric dt.mif -vector eig.mif');
    
    # compute pdd map from eig and fa maps
    os.system('mrcalc eig.mif fa.mif -mult pdd.mif');
    
    # compute response function from dwi data for CSD
    os.system('dwi2response dwi.mif rsp.txt -mask mask.mif -grad grads -nthreads 12');
    
    # estimate fiber orientation distribution from response and DWI, then track
    # fibers csd max number of fibers 1000000
    os.system('dwi2fod dwi.mif rsp.txt FOD.mif -mask mask.mif -grad grads -nthreads 12');
    
    # perform probabilistic tracography on FOD fields
    os.system('tckgen FOD.mif full_brain.tck -seed_image dwi.mif -mask mask.mif -number 1000000 -nthreads 12');
    