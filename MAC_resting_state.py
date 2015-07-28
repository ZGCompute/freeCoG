# MAC_resting_state.py
# Author: Zack Greenberg
# UCSF Neuroimaging Center
# Department of Neurology
# Last Editd May 1st, 2014

# This script contains the functions/modules called by run_RESTpipe.py to
# run the MAC pipeline for BOLD resting-state proessing.

import re
import shutil
import glob
import tempfile
import logging
from os import path
import subprocess
import os
import csv
import getpass
import nipype
import nipype.interfaces.fsl as fsl
import nipype.interfaces.spm as spm
from nipype.interfaces.spm import SliceTiming
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.matlab as mlab
from time import sleep
import signal
from datetime import datetime
from random import choice
from zipfile import ZipFile as zf


VERSION = "3.7.0";
PROCESS_NAME = "Resting State";

_log = logging.getLogger();

# Get fsl version and freesurfer version
os.environ['FSLVERSION'] = '5.0.5';
os.environ['FSVERSION'] = '5.1';

# where asl matlab scripts are kept
matlab_dir= '/usr/local/MATLAB/R2014a';
_matlab_scripts_path = os.path.normpath(os.path.join(matlab_dir));

# where spm8 is kept, on cloudi its /home/zgreenberg/Desktop/spm8
_spm_path = os.path.join('/usr/local/MATLAB/R2014a/toolbox/spm8');
mlab.MatlabCommand.set_default_paths(_spm_path);

def get_project_dir(group_name, project_name):

    #The path to the root folder for the project where Rest_Pipe is located
    os.chdir('/mnt/macdata/groups');
    root_dir= os.getcwd();
    project_dir = root_dir + '/' + group_name + '/' + project_name;
    
    # Throw exception if Rest-Pipe folder is not found
    if not path.exists(project_dir):
        raise Exception("Project folder, %s, does not exist! Please request that it be created by your system administrator" % project_dir)

    return project_dir
    

        
def get_case_dir(project_dir, accessionNum):

    # Function finds the path to the root folder for a given case
    # This is the user entered accession number
    a = accessionNum;

    # Find the folder in the Rest_Pipe dir containing the subject's accession number
    for fname in os.listdir(project_dir):
        if fname.startswith(a):
            patientFolder = fname;
    
    # Create a variable for the patient's folder directory, error if it doesn't exist
    patient_dir = project_dir + '/' + patientFolder;
    if not path.exists(patient_dir):
       raise Exception, "Cannot figure out the working dir for patient: %s Code: %s" % (project_name, a)
        
    return patient_dir
        
        
                
def find_asl_t2(case_dir):

   # Function searches for the AC-PC Aligned T2 image that will be used for Distortion Correction and Coregistration
   # Find the T2 nifti file    
   for fname in os.listdir(case_dir):
      if fname.startswith('T2_') and fname.endswith('.nii'):
         TFile = fname;
 
   # Go into the patients folder, make a dir for the ACPC aligned T2 and move it into that dr       
   os.chdir(case_dir)
   os.mkdir(case_dir + '/' + 'ACPC_Alignment')
   ACPC_Alignment= case_dir + '/' + 'ACPC_Alignment'
   shutil.move(TFile, ACPC_Alignment);
   ACPC_dir = case_dir + '/' + 'ACPC_Alignment' + '/' + TFile;

   return ACPC_dir;  


def find_rest_t1(case_dir):

   """Search for the AC-PC Aligned T1 image that will be used for segmentation and partial volume correction"""
   
   os.chdir(case_dir)
   os.mkdir('T1_Segmentation');
   T1_dir = case_dir + '/' + 'T1_Segmentation';     
   for fname in os.listdir(case_dir):

      if fname.startswith('MP-LAS-long') and fname.endswith('.nii'):
         T1File = fname;
         shutil.move(T1File, T1_dir);
         break;
      elif fname.startswith('MP-LAS-long') and fname.endswith('.zip'):
         T1File = fname;
         shutil.move(T1File, T1_dir);
         os.chdir(T1_dir);
         with zf(T1File) as zf_name:
            zf_name.extractall();
         os.system('dcm2nii *');
         break;
      elif fname.startswith('MP-LAS') and fname.endswith('.nii'):
         T1File = fname;
         shutil.move(T1File, T1_dir);
         break;
      elif fname.startswith('MP-LAS') and fname.endswith('.zip'):
         T1File = fname;
         shutil.move(T1File, T1_dir);
         os.chdir(T1_dir);
         with zf(T1File) as zf_name:
            zf_name.extractall();
         os.system('dcm2nii *');
         break;

   return T1_dir;   

def get_rest_rawz_dir(case_dir):


   # Funtion finds the zipped rest-raw folder and returns the directory name to unzip it	
   for fname in os.listdir(case_dir):
      if fname.startswith('rsfMRI') and fname.endswith('.zip'):
         rest_rawz_dir = case_dir + '/' + fname;

   return rest_rawz_dir


def unzip_rest_get_dir(rest_rawz_dir, case_dir):

   # Function unzips the raw asl dicom dir and returns new unzipped dir.
   with zf(rest_rawz_dir) as zf_dir:
      zf_dir.extractall('rsfMRI-raw')

   for fname in os.listdir(case_dir + '/' + 'rsfMRI-raw'):
      rest_sub_dir = fname;
   
   # Create a variable for the dicom directory
   dicomdir = case_dir + '/' + 'rsfMRI-raw' + '/' + rest_sub_dir;

   return dicomdir; 

def convert_dcm_tonii(dicom_dir):

   # Function Converts dcm EPIs to .nii  
   os.chdir(dicom_dir);
   os.mkdir('nii_all');
   for files in os.listdir(os.getcwd()):
      os.system('dcm2nii -a n -d n -e n -g n -i n -p n -f y -n y -v n %s' %(files));

   os.system('mv *.nii nii_all');
   nii_dir = dicom_dir + '/nii_all';

   return nii_dir;

def run_bet(nii_dir):

   # Function uses FSL Bet to skullstrip all the .nii EPIs   
   os.chdir(nii_dir);
   for files in os.listdir(os.getcwd()):
      if files.endswith('.nii'):
         btr = fsl.BET();
         btr.inputs.in_file = files;
	 btr.inputs.frac= 0.7;
	 print '<:::Extracting brain from EPI %s:::>' %(files)
	 res = btr.run();

def run_spm_realign(nii_dir):

   # Function realigns all EPIs to first aquisition in sequence 
   # Create a container for names of skull-stripped EPIs
   stripped_list =[];
   unzipped_stripped_list = [];
   
   # Move into tagegd dir, get a list of all skull-stripped files, and create sub dir for stripped images
   os.chdir(nii_dir);
   os.mkdir('skull_stripped');
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii.gz'):
         stripped_list.append(files);

   # Sort list
   stripped_list.sort();
   for fname in stripped_list:
      os.system('mv %s skull_stripped' %(fname));

   # Unzip the .gz files
   os.chdir('skull_stripped');
   for fname in stripped_list:
      new_fname = fname[:-3];
      os.system('gunzip %s %s' %(fname, new_fname));

   # Get final list of unzipped skull-stripped files
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii'):
         unzipped_stripped_list.append(files);

   # Run spm realign on tagged skull stripepd images
   unzipped_stripped_list.sort();   
   realign = spm.Realign();
   realign.inputs.in_files = unzipped_stripped_list;
   realign.inputs.register_to_mean = False;
   realign.inputs.paths = os.getcwd();
   print '<:::Realigning tagged epis to first aquisition in sequence:::>';
   realign.run();
   stripped_aligned_dir = os.getcwd();

   return stripped_aligned_dir;


def run_slice_time_correction(stripped_aligned_dir):

   # Function runs slice time correction for an interleived ascending EPI sequence
   # Navigate to dir with skull-stripped/realigned EPIs, create a list of the filenames to iter over
   os.chdir(stripped_aligned_dir);
   aligned_list = [];
   for files in os.listdir(os.getcwd()):
      if files.startswith('rIM'):
         aligned_list.append(files);

   # Initialize and run spm slicetiming correction
   aligned_list.sort();
   st = SliceTiming();
   st.inputs.in_files = aligned_list;
   st.inputs.num_slices = 36;
   st.inputs.time_repetition = 2.0;
   st.inputs.time_acquisition = 2. - (2./36.);
   st.inputs.slice_order = range(1,37,2) + range(2,38,2);
   st.inputs.ref_slice = 18;
   st.inputs.paths = os.getcwd();
   print "<:::Correcting EPIs for slice-timeing:::>";
   st.run();

   os.mkdir('slice_timed');
   os.system('mv ar* slice_timed');
   os.chdir('slice_timed');
   slice_time_dir = os.getcwd();

   return slice_time_dir;
     
   
def coregister_to_T1(slice_time_dir, T1_dir):

   # Function coregisters all EPIs in sequence to ac-pc aligned T1
   # Get the name of the ac-pc aligned t1
   os.chdir(T1_dir);
   for files in os.listdir(T1_dir):
      T1name = files;

   # Go to the dir with the slice time corrected EPI's and iter over 
   # filenames to coregister each to the ac-pc aligned T1.
   os.chdir(slice_time_dir);
   for files in os.listdir(slice_time_dir):
      if files.startswith("ar"):
         fname = files;
         coreg = spm.Coregister();
         coreg.inputs.target = slice_time_dir + '/' + fname;
         coreg.inputs.source = T1_dir + '/' + T1name;
         print "<:::Coregistering EPI %s to AC-PC aligned T1:::>" %(fname);
         coreg.run();

def normalize_T1_toMNI152(T1_dir,slice_time_dir):

   # Function normalizes an AC-PC aligned T1 into MNI 152 space
   # Get the name of the ac-pc aligned T1
   os.chdir(T1_dir);
   for files in os.listdir(T1_dir):
      if files.startswith('r'):
         T1name = files;
  
   # Get list of EPIs to apply MNI normalization transform result from T1
   EPI_list = [];
   os.chdir(slice_time_dir);
   for files in os.listdir(slice_time_dir):
      if files.startswith('ar'):
         EPI_list.append(files);

   # Initialize SPM normalization
   EPI_list.sort();
   mni_template = matlab_dir + '/toolbox/spm8/canonical/avg152T1.nii';
   os.system('cp %s %s' %(mni_template,slice_time_dir)); # copy spm mni template to cwd
   norm = spm.Normalize();
   norm.inputs.source = T1_dir + '/' + T1name;
   norm.inputs.template = mni_template;
   norm.inputs.affine_regularization_type = 'mni';
   EPI_list.append(norm.inputs.source); # apply normalization to T1 and all EPIs
   norm.inputs.apply_to_files = EPI_list;
   print "<:::Normalizing T1 and EPIs to MNI 152:::>";
   norm.run();

   # Make a new dir for MNI warped EPIs 
   MNI_dir = os.mkdir('MNI_Normalized');
   os.system('cp w* %s' %(MNI_dir));

   return MNI_dir;

 
def extract_EPI_timeseries(MNI_dir, T1_dir):

   # Function merges all MNI normalized EPIs into single 4d vol (FSL)
   # extracts a time series from MNI ROIs
   os.chdir(MNI_dir);
   os.system('fslmerge -t rsfMRI.nii w*');
   
   
