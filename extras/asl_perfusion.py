# asl_perfusion.py
# Author: Zack Greenberg, addapted from Derek Flenniken
# UCSF Neuroimaging Center
# Department of Neurology
# Last Editd May 1st, 2014

# This script contains the functions/modules called by run_ASLpipe.py to
# run the MAC pipeline for pulsed ASL proessing.

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
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.matlab as mlab
from time import sleep
import signal
from datetime import datetime
from random import choice
from zipfile import ZipFile as zf


VERSION = "3.7.0";
PROCESS_NAME = "ASL Perfusion";

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

    #The path to the root folder for the project where ASL-pipe is located
    os.chdir('/mnt/macdata/groups');
    root_dir= os.getcwd();
    project_dir = root_dir + '/' + group_name + '/' + project_name;
    
    # Throw exception if ASL-pipe folder is not found
    if not path.exists(project_dir):
        raise Exception("Project folder, %s, does not exist! Please request that it be created by your system administrator" % project_dir)

    return project_dir
    

        
def get_case_dir(project_dir, accessionNum):

    # Function finds the path to the root folder for a given case
    # This is the user entered accession number
    a = accessionNum;

    # Find the folder in the ASL_pipe dir containing the subject's accession number
    for fname in os.listdir(project_dir):
        if fname.startswith(a):
            patientFiles = fname;
    
    # Create a variable for the patient's folder directory, error if it doesn't exist
    patient_dir = project_dir + '/' + patientFiles;
    if not path.exists(patient_dir):
        raise Exception, "Cannot figure out the working dir for patient: %s Code: %s" % (project_name, a)
        
    return patient_dir
        
        
                
def find_asl_t2(case_dir):

   # Function searches for the AC-PC Aligned T2 image that will be used for Distortion Correction and Coregistration
   # Find the T2 nifti file    
   for fname in os.listdir(case_dir):
      if fname.startswith('T2_') and fname.endswith('.nii'):
         TFile = fname;
 
   # Go into the patients fodler, make a dir for the ACPC aligned T2 and move it into that dr       
   os.chdir(case_dir)
   os.mkdir(case_dir + '/' + 'ACPC_Alignment')
   ACPC_Alignment= case_dir + '/' + 'ACPC_Alignment'
   shutil.move(TFile, ACPC_Alignment);
   ACPC_dir = case_dir + '/' + 'ACPC_Alignment' + '/' + TFile;

   return ACPC_dir;  


def get_asl_rawz_dir(case_dir):


   # Funtion finds the zipped ASL-raw folder and returns the directory name to unzip it	
   for fname in os.listdir(case_dir):
      if fname.startswith('ASL-raw') and fname.endswith('.zip'):
         asl_rawz_dir = case_dir + '/' + fname;

   return asl_rawz_dir

def unzip_asl_get_dir(asl_rawz_dir, case_dir):

   # Function unzips the raw asl dicom dir and returns new unzipped dir.
   with zf(asl_rawz_dir) as zf_dir:
      zf_dir.extractall('ASL-raw')

   for fname in os.listdir(case_dir + '/' + 'ASL-raw'):
      asl_sub_dir = fname;
   
   # Create a variable for the dicom directory
   dicomdir = case_dir + '/' + 'ASL-raw' + '/' + asl_sub_dir;

   return dicomdir; 

def run_PWI_procAll(dicom_dir):

   # Make a directory for all nifitis (tagged and untagged), 
   # convert all the dcms to nii, move them to the nii_all dir,
   # skull-strip them, realign them to first EPI, so we can use them
   # in fsl's asl-subtract to get perfusion weighted image (PWI avg).
   os.chdir(dicom_dir)
   os.mkdir('nii_all');
   for files in os.listdir(os.getcwd()):
      os.system('dcm2nii -a n -d n -e n -g n -i n -p n -f y -n y -v n %s' %(files));

   # move into the nii_all dir and skull-strip/realign
   os.system('mv *.nii nii_all');
   os.chdir('nii_all');
   
   # Run FSL BET()
   for files in os.listdir(os.getcwd()):
      if files.endswith('.nii'):
         btr = fsl.BET();
         btr.inputs.in_file = files;
	 btr.inputs.frac= 0.6;
         btr.inputs.robust = True;
	 print 'Extracting brain from epi %s......' %(files)
	 res = btr.run();
 
   # Make skull-stripped_realigned dir
   os.mkdir('realigned_stripped');
   os.system('gunzip *brain.nii.gz');
   
   # Get list of all stripped EPIs in dir
   #os.system('rm mean*');
   stripped_list =[];
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii'):
         stripped_list.append(files);
  
   # Sort the list to be realigned and remove m0 from the end
   stripped_list.sort();
   m0 = stripped_list.pop();
   if m0 != "m0_brain.nii":
      print "Error: EPI %s was incorectly excluded from realignment"

   # Run spm_realign to realign EPIs to first in sequence   
   realign = spm.Realign();
   realign.inputs.in_files = stripped_list;
   realign.inputs.register_to_mean = False;
   realign.inputs.paths = os.getcwd();
   print 'Realigning tagged and untagged epis to m0...';
   realign.run();

   # Move the stripped and realigned EPIs to seperate dir
   os.system('mv r* realigned_stripped');
   os.chdir('realigned_stripped');
  
   # Rename first EPI aquisition m0 (non PWI image)
   # and create 4D nii to run asl_subtract on: 
   # (subtracts even-untagged from odd-tagged volumes).
   os.system('rm *001_brain.nii'); # remove first volume (m0) before performing subtraction
   os.system('rm *.txt'); # Remove text file before combining all nii into 4d
   os.system('fslmerge -t asl.nii r*');
   os.system('asl_file --data=asl.nii --ntis=1 --iaf=tc --diff --out=diffdata --mean=diffdata_mean');
   os.system('gunzip *.nii.gz');
   all_aligned_dir = os.getcwd();
 
   return all_aligned_dir;

   
def CBFscale_PWI_data(all_aligned_dir):

   # Envoke matlab to calculate scaled-CBF from PWI data
   mlc = mlab.MatlabCommand();
   cmd = "cd('%s');raw_pwi = spm_vol('diffdata_mean.nii');scaled_pwi = raw_pwi;scaled_pwi.fname = 'CBF_Scaled_PWI.nii';scaled_pwi.descript = 'Scaled from the PWI Image';pwi_data = spm_read_vols(raw_pwi);Lamda = 0.9000;Alpha = 0.9500;Tau = 22.50;R1A = (1684)^-1;PER100G = 100;SEC_PER_MIN = 60;MSEC_PER_SEC = 1000;TI1 = 700;TI2 = 1800;PWI_scale = zeros(size(pwi_data));sliceNumbers = (1:size(pwi_data, 3))';Constant = Lamda / (2 * Alpha * TI1) * (PER100G * SEC_PER_MIN * MSEC_PER_SEC);Slice_based_const = exp(R1A * (TI2 + (sliceNumbers - 1) * Tau));Numerator = pwi_data;for n =1:size(sliceNumbers);   PWI_scale(:,:,n) = Constant * Slice_based_const(n) * Numerator(:,:,n);end; spm_write_vol(scaled_pwi, PWI_scale);" %(all_aligned_dir);

   mlc.inputs.script = cmd;
   out = mlc.run();

                
def sort_tagged_untagged(dicom_dir):

   # Function finds even numbered(untagged), odd numbered(tagged), and m0 EPIs and sort/rename them.  
   os.chdir(dicom_dir);
   os.mkdir('tagged');
   os.mkdir('untagged');

   # Place even numbered aquistions in untagged folder, and odd aquisitions in tagged folder
   for files in os.listdir(dicom_dir):

      if files == 'tagged' or files == 'untagged':
         pass; # skipp dir names

      elif float(files[9:12])%2 == 0 and files[9:12] != '001':
         dir = os.getcwd();new_fileName = dir + '/' + 'untagged/untagged_' + files[9:12];
         os.system('cp %s %s' %(files,new_fileName));

      elif float(files[9:12])%2 != 0 and files[9:12] != '001':
         dir = os.getcwd();new_fileName = dir + '/' + 'tagged/tagged_' + files[9:12];
         os.system('cp %s %s' %(files,new_fileName));
	
      elif files[9:12] == '001':
	 dir = os.getcwd();
         new_fileName = dir + '/' + 'm0'; # create m0 from first non-perfusion weighted epi
         os.system('cp %s %s' %(files,new_fileName));

   # Store variables for tagged and untagged dirs, move a copy of m0 to each
   dir = os.getcwd();
   tagged_dir = dir + '/' + 'tagged/';
   untagged_dir = dir + '/' + 'untagged/';
   os.system('sudo cp m0 %s' %(tagged_dir));
   os.system('sudo cp m0 %s' %(untagged_dir));


   return tagged_dir,untagged_dir

def convert_tagged_untagged_tonii(tagged_dir,untagged_dir):

   # Function Converts tagged and untagged EPIs to .nii  
   os.chdir(tagged_dir);
   for files in os.listdir(os.getcwd()):
      if files.startswith('tagged_'):
         os.system('dcm2nii -a n -d n -e n -g n -i n -p n -f y -n y -v n %s' %(files));
      elif files.startswith('m0'):
         os.system('dcm2nii -a n -d n -e n -g n -i n -p n -f y -n y -v n %s' %(files));

   os.chdir(untagged_dir);
   for files in os.listdir(os.getcwd()):
      if files.startswith('untagged_'):
         os.system('dcm2nii -a n -d n -e n -g n -i n -p n -f y -n y -v n %s' %(files));
      elif files.startswith('m0'):
         os.system('dcm2nii -a n -d n -e n -g n -i n -p n -f y -n y -v n %s' %(files));



def run_bet(tagged_dir,untagged_dir):

   # Function uses FSL Bet to skullstrip all the .nii EPIs   
   os.chdir(tagged_dir);
   for files in os.listdir(os.getcwd()):
      if files.endswith('.nii'):
         btr = fsl.BET();
         btr.inputs.in_file = files;
	 btr.inputs.frac= 0.7;
	 print 'Extracting brain from epi %s......' %(files)
	 res = btr.run();

   os.chdir(untagged_dir);
   for files in os.listdir(os.getcwd()):
      if files.endswith('.nii'):
         btr = fsl.BET();
         btr.inputs.in_file = files;
	 btr.inputs.frac= 0.7;
	 print 'Extracting brain from epi %s......' %(files)
	 res = btr.run();


def run_spm_realign(tagged_dir, untagged_dir):

   # Function realigns all the tagged/untagged EPIs to first run in sequence (i.e., m0).
   # Create containers for tagged/stripped, and untagged/stripped to feed to realign function
   tagged_list =[];
   tagged_stripped_list = [];
   untagged_list =[];
   untagged_stripped_list = [];

   # Move into tagegd dir, get a list of all skull-stripped files, and create sub dir for stripped images
   os.chdir(tagged_dir);
   os.mkdir('skull_stripped');
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii.gz'):
         tagged_list.append(files);

   # Sort so m0 is first in the list (aligns to first epi)
   tagged_list.sort();
   for fname in tagged_list:
      os.system('mv %s skull_stripped' %(fname));

   # Unzip the .gz files
   os.chdir('skull_stripped');
   for fname in tagged_list:
      new_fname = fname[:-3];
      os.system('gunzip %s %s' %(fname, new_fname));

   # Get final list of unzipped skull-stripped files
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii'):
         tagged_stripped_list.append(files);

   # Run spm realign on tagged skull stripepd images
   tagged_stripped_list.sort();   
   realign = spm.Realign();
   realign.inputs.in_files = tagged_stripped_list;
   realign.inputs.register_to_mean = False;
   realign.inputs.paths = os.getcwd();
   print 'Realigning tagged epis to m0...';
   realign.run();
   tagged_aligned_dir = os.getcwd();

   # Move into untagegd dir, get a list of all skull-stripped files, and create sub dir for stripped images
   os.chdir(untagged_dir);
   os.mkdir('skull_stripped');
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii.gz'):
         untagged_list.append(files);

   # Sort so m0 is first in the list (aligns to first epi)
   untagged_list.sort();
   for fname in untagged_list:
      os.system('mv %s skull_stripped' %(fname));

   # Unzip the .gz files
   os.chdir('skull_stripped');
   for fname in untagged_list:
      new_fname = fname[:-3];
      os.system('gunzip %s %s' %(fname, new_fname));

   # Get final list of unzipped skull-stripped files
   for files in os.listdir(os.getcwd()):
      if files.endswith('brain.nii'):
         untagged_stripped_list.append(files);

   # Run spm realign on untagged skull stripepd images
   untagged_stripped_list.sort();
   realign = spm.Realign();
   realign.inputs.in_files = untagged_stripped_list;
   realign.inputs.register_to_mean = False;
   print 'Realigning untagged epis to m0...'; 
   realign.run();
   untagged_aligned_dir = os.getcwd();

   return tagged_aligned_dir, untagged_aligned_dir;



def Get_perfusion_calc(case_dir, tagged_aligned_dir, untagged_aligned_dir):
   
   # Function sums and avgs skull stripped/aligned EPIs for tagged and untagged aquisitions """

   raw_perfusion_dir = os.mkdir(case_dir + '/Raw_Perfusion');
   aligned_list =[];
   os.chdir(tagged_aligned_dir);
   for files in os.listdir(os.getcwd()):
      if files.startswith('rtagged'):
         aligned_list.append(files);

   aligned_list.sort();
   maths = fsl.ImageMaths(in_file = aligned_list[0], op_string = '-add %s' %(aligned_list[1]), out_file = 'tagged_sum.nii.gz')
   maths.run();
   for fname in aligned_list[2:]:
      print 'Summing tagged EPI %s' %(fname)
      maths = fsl.ImageMaths(in_file = fname, op_string = '-add %s' %('tagged_sum.nii.gz'), out_file = 'tagged_sum.nii.gz')
      maths.run();

   denom = len(aligned_list);
   maths = fsl.ImageMaths(in_file = 'tagged_sum.nii.gz', op_string = '-div %s' %(denom), out_file = 'tagged_avg.nii.gz')
   maths.run();
   os.system('mv tagged_avg.nii.gz %s' %(raw_perfusion_dir));
   aligned_list =[];
   os.chdir(untagged_aligned_dir);
   for files in os.listdir(os.getcwd()):
      if files.startswith('runtagged'):
         aligned_list.append(files);

   aligned_list.sort();
   maths = fsl.ImageMaths(in_file = aligned_list[0], op_string = '-add %s' %(aligned_list[1]), out_file = 'untagged_sum.nii.gz')
   maths.run();
   for fname in aligned_list[2:]:
      print 'Summing tagged EPI %s' %(fname)
      maths = fsl.ImageMaths(in_file = fname, op_string = '-add %s' %('untagged_sum.nii.gz'), out_file = 'untagged_sum.nii.gz')
      maths.run();

   denom = len(aligned_list);
   maths = fsl.ImageMaths(in_file = 'untagged_sum.nii.gz', op_string = '-div %s' %(denom), out_file = 'untagged_avg.nii.gz')
   maths.run();
   os.system('mv untagged_avg.nii.gz %s' %(raw_perfusion_dir));

   os.chdir(raw_perfusion_dir);
   maths = fsl.ImageMaths(in_file = 'tagged_avg.nii.gz', op_string = '-sub %s' %('untagged_avg.nii.gz'), out_file = 'mean_perfusion_raw.nii.gz')
   maths.run();

   return raw_perfusion_dir;

   

def find_asl_t1(case_dir):

   """Search for the AC-PC Aligned T1 image that will be used for segmentation and partial volume correction"""
   
   os.chdir(case_dir)
   os.mkdir('PVE_Segmentation');
   PVE_Segmentation_dir = case_dir + '/' + 'PVE_Segmentation';     
   for fname in os.listdir(case_dir):

      if fname.startswith('MP-LAS-long') and fname.endswith('.nii'):
         T1File = fname;
         shutil.move(T1File, PVE_Segmentation_dir);
         break;
      elif fname.startswith('MP-LAS-long') and fname.endswith('.zip'):
         T1File = fname;
         shutil.move(T1File, PVE_Segmentation_dir);
         os.chdir(PVE_Segmentation);
         with zf(T1File) as zf_name:
            zf_name.extractall();
         os.system('dcm2nii *');
         break;
      elif fname.startswith('MP-LAS') and fname.endswith('.nii'):
         T1File = fname;
         shutil.move(T1File, PVE_Segmentation_dir);
         break;
      elif fname.startswith('MP-LAS') and fname.endswith('.zip'):
         T1File = fname;
         shutil.move(T1File, PVE_Segmentation_dir);
         os.chdir(PVE_Segmentation);
         with zf(T1File) as zf_name:
            zf_name.extractall();
         os.system('dcm2nii *');
         break;

   return PVE_Segmentation_dir;   



def run_spm_segmentT1(PVE_segmentation_dir):

   # Go into dir with acpc_aligned t1 and find the t1 filename
   os.chdir(PVE_segmentation_dir);
   for files in os.listdir(os.getcwd()):
      if files.endswith(".nii"):
         T1 = files;
   
   # Run Spm_NewSegment on the T1 to get Gm,wm,ventricles
   print "Segmenting Grey matter and White matter from T1 structural...";
   seg = spm.NewSegment();
   seg.inputs.channel_files = T1;
   seg.inputs.channel_info = (0.0001, 60, (True, True))
   seg.run();
