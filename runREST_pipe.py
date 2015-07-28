# runREST_pipe.py

# Description: This script imports a set of python modules to be used
# for processing fMRI BOLD connecitivity data from the UCSF Neuroimaging Center and 
# the UCSF Memory and Aging Center. The script should be run on the cloud at:
# /mnt/macdata/groups/Rest_Pipe/. If you drop your data in that directory
# and call the fallowing commands, this script will run: 

# 1.) Dicom to niftii conversion (dcm2nii)
# 2.) Skull stripping (FSL BET)
# 3.) Realigning the EPIs to the first aquisition in EPI sequence.(SPM)
# 4.) Register Ac-Pc aligned T2 weighted image to EPIs (SPM)
# 5.) Run distortion correction using T2 image (MATLAB, ITK)
# 6.) Register your EPIs to an Ac-Pc aligned T1 (SPM)
# 7.) Run Partial Volume Correction (MATLAB)
# 8.) Normalize EPIs to MNI 152 space (SPM)
# 9.) Extract a BOLD time series from MNI defined ROIs (FSL)
# 10.) Filter the time-series at 0.08 - 0.1Hz (Matlab)
# 11.) Run pairwise partial-corellations between all ROIs.(Matlab)
# 12.) Run a GLM for specific seed ROIs (DMN, Motor, Language), include motion, GMWM, Ventricles as regressors (Matlab)
# 13.) Correct for multiple comparisons using a surrogate simulation (MATLAB)

# Dependencies: Matlab, Python 2.7 >, Nipype, SPM, FSL, Freesurfer, ITK

# Author: Zack Greenberg
# Department of Neurology and Physiology
# University of California, San Francisco
# Last Edited: May 15th, 2014

#resting-state repository
from MAC_resting_state import * 

# Folders to naviagate to on /mnt/macdata/groups/
group_name = 'imaging_core';
project_name='Rest_Pipe';
accessionNum= raw_input("Welcome to MAC Pipe! Please enter your subject's accession number:  "); 

# Define root and data dir variables
project_dir =get_project_dir(group_name,project_name); 
case_dir=get_case_dir(project_dir,accessionNum);
ACPC_dir = find_asl_t2(case_dir); # T2 used for coreg and distortion correction
T1_dir = find_rest_t1(case_dir); # T1 for segmenting gmwm for partial volume correction 
rest_rawz_dir = get_rest_rawz_dir(case_dir);

# 1.) Unzip raw dicom EPIs and return the dir
dicom_dir = unzip_asl_get_dir(asl_rawz_dir, case_dir);

# 2.) dicom to niftii conversion
nii_dir = convert_dcm_tonii(dicom_dir);

# 3.) Skull stripping using FSL bet
run_bet(nii_dir);

# 4.) Realigning the EPIs to the first EPI using spm_realign (Motion-correction). 
stripped_aligned_dir = run_spm_realign(nii_dir);

# 5.) Run spm slice-time correction on skulstripped/aligned EPIs. 
slice_time_dir = run_slice_time_correction(stripped_aligned_dir);

# 6.) Run spm coregistration to bring EPIs into T1 space
coreg_dir = coregister_to_T1(slice_time_dir, T1_dir);

# 7.) Normalize EPIs to MNI 152 space (SPM)
MNI_dir = normalize_T1_toMNI152(T1_dir,slice_time_dir);
