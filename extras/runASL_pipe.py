# runASL_pipe.py

# Description: This script imports a set of python modules to be used
# for processing pulsed ASL perfusion data from the UCSF Neuroimaging Center and 
# the UCSF Memory and Aging Center. The script should be run on the cloud at:
# /mnt/macdata/groups/ASL_pipe/. If you drop your data in that directory
# and call the fallowing commands, this script will run: 

# 1.) dicom to niftii conversion (dcm2nii)
# 2.) Sorting of tagged (perfusion weighted) and untagged EPIs
# 3.) Skull stripping (FSL BET)
# 4.) Realigning the EPIs to the non-perfusion weighted m0 (i.e., first aquisition in EPI sequence) (SPM)
# 5.) Calculate ASL perfusion maps via subtraction of mean images from tagged and untagged (FSL)
# 6.) Register Ac-Pc aligned T2 weighted image to EPIs (SPM)
# 7.) Run distortion correction of perfusion maps using T2 image (MATLAB)
# 8.) Register ASL maps to T2, and T2 to T1, then combine the affine transforms (SPM)
# 9.) Run Partial Volume Correction (MATLAB)
# 10.) Calculate and Normalize CBF maps (CBF/m0) (MATLAB)
# 11.) Extract mean CBF from freesurfer defined ROIs (FREESURFER)

# Author: Zack Greenberg, addapted from Derek Flenelen
# Department of Neurology and Physiology
# University of California, San Francisco
# Last Edited: April 16th, 2014

#ASL repository
from MAC_asl_perfusion import * 

# Folders to naviagate to on /mnt/macdata/groups/
group_name = 'imaging_core';
project_name='ASL_pipe';
accessionNum= raw_input("Please enter subject's accession number:  "); #'PPG0246-1' for test; 

# Define root and data dir variables
project_dir =get_project_dir(group_name,project_name); 
case_dir=get_case_dir(project_dir,accessionNum);
ACPC_dir = find_asl_t2(case_dir); # T2 used for coreg and distortion correction
PVE_Segmentation_dir = find_asl_t1(case_dir); # T1 for segmenting gmwm for partial volume correction 
asl_rawz_dir = get_asl_rawz_dir(case_dir);

# 1.) Sorting of tagged (perfusion weighted) and untagged EPIs
dicom_dir = unzip_asl_get_dir(asl_rawz_dir, case_dir);
[tagged_dir,untagged_dir] = sort_tagged_untagged(dicom_dir);

# 2.) Run dcm2nii conversion, skull-stripping(FSL BET), realignment(SPM), and PWI calc (FSL)
# This runs on all dicoms/niis without sorting by tagged/untagged first.
all_aligned_dir = run_PWI_procAll(dicom_dir);

# 3.) Scale PWI for time lag and compute CBF
CBFscale_PWI_data(all_aligned_dir);

# 4.) dicom to niftii conversion
convert_tagged_untagged_tonii(tagged_dir,untagged_dir);

# 5.) Skull stripping using FSL bet
run_bet(tagged_dir,untagged_dir);

# 6.) Realigning the EPIs to the non-perfusion weighted m0 using spm_realign. 
[tagged_aligned_dir, untagged_aligned_dir] = run_spm_realign(tagged_dir, untagged_dir);

# 7.) Calculate ASL perfusion maps via subtraction of mean images from tagged and untagged aquisitions
Get_perfusion_calc(case_dir, tagged_aligned_dir, untagged_aligned_dir);

# 8.) Run spm_NewSegment to get GMWM and ventricles for Partial volume correction
run_spm_segmentT1(PVE_segmentation_dir);
