# img_pipe.py

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: July 21, 2015

# This file contains the chang lab imaging pipeline (freeCoG)
# as one importable python class for running a patients
# brain surface reconstruction and electrode localization/labeling
# as well as dwi and tractography for reconstructing 
# major white matter bundles and full-brain connectome

import os
import numpy as np
import scipy
import scipy.io
import nibabel as nib
import nipy
from scipy.spatial.distance import cdist
import collections
import time
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mahotas
import subprocess
import itertools
import shlex
from dipy.segment.mask import median_otsu
from nipy.algorithms.kernel_smooth import LinearFilter
import mayavi
from mayavi import mlab
from nipype.interfaces import matlab as matlab
from nipype.interfaces import spm
from nipype.interfaces import fsl
from nipype.interfaces.fsl import FLIRT



class freeCoG:

   # blank var for patient name
   subj = '';  
   
   # blank var for patient hem of implantation
   hem='';

   # fs subjects dir on server where you run img_pipe
   subj_dir = "/path/to/subjects/";
   
   # IP of home computer to requester data from
   home = "";
   
   # directory where imaging data is on your home machine
   # (not the server your running img_pipe on
   home_data_dir = '/path/to/subjects/';
   
   # set python and bash "areal" scripts dir on current server
   bash_dir = '/path/to/bash';
   python_dir = '/path/to/python';
   AFQ_dir = '/path/to/AFQ';
   SPM_dir = '/path/to/spm12';
   VISTA_dir = '/path/to/vistasoft-master';
   Hough_dir = '';
   project_elecs_dir = '/path/to/electrode_placement/';  
 
   # set dir for CT img data
   CT_dir = subj_dir + subj + '/CT/';
   
   # set dir for elecs coordinates
   elecs_dir = subj_dir + subj + '/elecs/';
   
   # set dir for DWI raw data and processed/tractography
   DWI_dir = '';
   DWI_bval =2000;
   
   # dir for fs look up table
   fsLUT_dir = '/usr/local/freesurfer/';
   fs_dir = subj_dir + '/' + subj + '/';


   # intialize class img_pipe constructor
   def __init__(self, subj_dir,subj):
      self.subj = subj;
      self.subj_dir = subj_dir;
      self.CT_dir = self.subj_dir + '/' + self.subj + '/CT/';
      self.elecs_dir = subj_dir + '/' + self.subj + '/elecs/';
      self.DWI_dir = subj_dir + '/' + self.subj + '/DTI_HARDI_55_2000/';
      self.fs_data_dir = subj_dir + '/' + self.subj + '/fs_data/';



   # method for running surface reconstruction using freesurfer
   def prep_recon(self):

	# nav to dir with subject fs folders 
	# and make new folder for patient
	os.chdir(self.subj_dir);
	os.mkdir(self.subj);
	os.chdir(self.subj);

	# creat mri and orig folders
	os.mkdir("mri");
	os.chdir("mri");
	os.mkdir("orig");
	os.chdir("orig");

	# grab subjects T1 from home machine
	os.system("scp %s:%s/%s/acpc/T1.nii ." %(self.home, self.home_data_dir, self.subj));

	# convert T1 to freesurfer 001.mgz format
	os.system("mri_convert T1.nii 001.mgz"); 

	# start tmux
	os.system("tmux");


   # function to run freesurfer recon-all
   def get_recon(self):
       
       # run recon-all for target subject
       os.system('recon-all -subjid %s -all -3T -openmp 12 -use-gpu' %(self.subj));
       
       

   # member function for converting freesurfer output .asc mesh files to .mat 
   # vertex and triangle coord lists
   def convert_fs2mlab(self,lh,rh):
       
      '''Member function for converting freesurfer output .asc mesh files to .mat
         vertex and triangle coord lists. lh and rh args are lh.pial and rh.pial
         freesurfer mesh file (full-path) respectively'''
    
    
      # use freesurfer mris_convert to get ascii pial surface for lh
      lh_ascii = lh + '.asc';
      os.system('mris_convert %s %s' %(lh, lh_ascii));
    
      # use freesurfer mris_convert to get ascii pial surface for rh
      rh_ascii = rh + '.asc';
      os.system('mris_convert %s %s' %(rh, rh_ascii));
    
      # clean up ascii file and extract matrix dimensions from header
      lh = open(lh_ascii, 'r');
      lh_mat = lh.readlines();
      lh.close();
      lh_mat.pop(0); # get rid of comments in header
      lh_mat = [item.strip('\n') for item in lh_mat];# get rid of new line char
    
      # extraxt inds for vert and tri    
      lh_inds = lh_mat.pop(0); 
      lh_inds = lh_inds.split(' '); # seperate inds into two strings
      lh_inds = [int(i) for i in lh_inds];# convert string to ints
    

      # get rows for vertices only, strip 0 column, and split into seperate strings
      lh_vert = [item.strip(' 0') for item in lh_mat[:lh_inds[0]]];
      lh_vert = [item.split('  ') for item in lh_vert]; # seperate strings
    
      # create containers for each column in vert matrix
      x = [];
      y = [];
      z = [];
    
      # fill the containers with float values of the strings in each column
      for i in lh_vert:
          x.append(float(i[0]));
          y.append(float(i[1]));
          z.append(float(i[2]));
        
      # convert to scipy mat
      x = scipy.mat(x);
      y = scipy.mat(y);
      z = scipy.mat(z);
    
      # concat columns to one n x 3 matrix
      x = x.transpose();
      y = y.transpose();
      z = z.transpose();
      lh_vert = scipy.column_stack((x,y,z));
      scipy.io.savemat('lh_vert.mat', {'vert':lh_vert});# save vert matrix
    
      # get rows for triangles only, strip 0 column, and split into seperate strings
      lh_tri = [item[:-2] for item in lh_mat[lh_inds[0] +1:]];
      lh_tri = [item.split(' ') for item in lh_tri]; # seperate strings  
    
      # create containers for each column in vert matrix
      x = [];
      y = [];
      z = [];
    
      # fill the containers with float values of the strings in each column
      for i in lh_tri:
         x.append(int(i[0]));
         y.append(int(i[1]));
         z.append(int(i[2]));
        
      # convert to scipy mat
      x = scipy.mat(x);
      y = scipy.mat(y);
      z = scipy.mat(z);
    
      # concat columns to one n x 3 matrix
      x = x.transpose();
      y = y.transpose();
      z = z.transpose();
      lh_tri = scipy.column_stack((x,y,z));
      lh_tri = lh_tri + 1;
      scipy.io.savemat('lh_tri.mat', {'tri':lh_tri});# save tri matrix  

      # convert inds to scipy mat
      lh_inds = scipy.mat(lh_inds);
      scipy.io.savemat('lh_inds.mat', {'inds':lh_inds});# save inds
    
      # repeat for rh
      rh = open(rh_ascii, 'r');
      rh_mat = rh.readlines();
      rh.close();
      rh_mat.pop(0); # get rid of comments in header
      rh_mat = [item.strip('\n') for item in rh_mat];# get rid of new line char
    
      # extraxt inds for vert and tri    
      rh_inds = rh_mat.pop(0); 
      rh_inds = rh_inds.split(' '); # seperate inds into two strings
      rh_inds = [int(i) for i in rh_inds];# convert string to ints
    

      # get rows for vertices only, strip 0 column, and split into seperate strings
      rh_vert = [item.strip(' 0') for item in rh_mat[:rh_inds[0]]];
      rh_vert = [item.split('  ') for item in rh_vert]; # seperate strings
    
      # create containers for each column in vert matrix
      x = [];
      y = [];
      z = [];
    
      # fill the containers with float values of the strings in each column
      for i in rh_vert:
         x.append(float(i[0]));
         y.append(float(i[1]));
         z.append(float(i[2]));
        
      # convert to scipy mat
      x = scipy.mat(x);
      y = scipy.mat(y);
      z = scipy.mat(z);
    
      # concat columns to one n x 3 matrix
      x = x.transpose();
      y = y.transpose();
      z = z.transpose();
      rh_vert = scipy.column_stack((x,y,z));
      scipy.io.savemat('rh_vert.mat', {'vert':rh_vert});# save vert matrix
    
      # get rows for triangles only, strip 0 column, and split into seperate strings
      rh_tri = [item[:-2] for item in rh_mat[rh_inds[0] +1:]];
      rh_tri = [item.split(' ') for item in rh_tri]; # seperate strings  
    
      # create containers for each column in vert matrix
      x = [];
      y = [];
      z = [];
    
      # fill the containers with float values of the strings in each column
      for i in rh_tri:
         x.append(int(i[0]));
         y.append(int(i[1]));
         z.append(int(i[2]));
        
      # convert to scipy mat
      x = scipy.mat(x);
      y = scipy.mat(y);
      z = scipy.mat(z);
    
      # concat columns to one n x 3 matrix
      x = x.transpose();
      y = y.transpose();
      z = z.transpose();
      rh_tri = scipy.column_stack((x,y,z));
      rh_tri = rh_tri + 1;
      scipy.io.savemat('rh_tri.mat', {'tri':rh_tri});# save tri matrix  

      # convert inds to scipy mat
      rh_inds = scipy.mat(rh_inds);
      scipy.io.savemat('rh_inds.mat', {'inds':rh_inds});# save inds



   # method for obtaining .mat files for vertex and triangle
   # coords of all subcortical freesurfer segmented meshes; needs self.subcortFs2mlab
   def get_subcort(self):

      '''Method for obtaining .mat files for vertex and triangle
         coords of all subcortical freesurfer segmented meshes'''
    
      # set ascii dir name
      subjAscii_dir = ('%s/%s/ascii/' %(self.subj_dir,self.subj));
    
      # tessellate all subjects freesurfer subcortical segmentations
      os.system('./aseg2srf.sh -s "%s" -l "4 5 10 11 12 13 17 18 26 \
                 28 43 44  49 50 51 52 53 54 58 60 14 15 16"' %(self.subj));
    
      # get list of all .srf files and change fname to .asc
      srf_list = [fname for fname in os.listdir(subjAscii_dir)];
      asc_list = [fname.replace('.srf', '.asc') for fname in srf_list];
      asc_list.sort();
      for fname in srf_list:
         new_fname = fname.replace('.srf', '.asc');
         os.system("mv %s %s" %(fname, new_fname));
    
      # convert all ascii subcortical meshes to matlab vert, tri coords
      subcort_list =['aseg_058.asc','aseg_054.asc', 'aseg_050.asc',\
                      'aseg_052.asc','aseg_053.asc', 'aseg_051.asc','aseg_049.asc',\
                       'aseg_043.asc', 'aseg_044.asc', 'aseg_060.asc', 'aseg_004.asc',\
                       'aseg_005.asc', 'aseg_010.asc', 'aseg_011.asc', 'aseg_012.asc',\
                       'aseg_013.asc', 'aseg_017.asc', 'aseg_018.asc', 'aseg_026.asc', \
                       'aseg_028.asc', 'aseg_014.asc', 'aseg_015.asc', 'aseg_016.asc'];
                       
      nuc_list = ['rAcumb', 'rAmgd', 'rCaud', 'rGP', 'rHipp','rPut', 'rThal',\
                   'rLatVent', 'rInfLatVent', 'rVentDienceph', 'lLatVent', 'lInfLatVent', \
                    'lThal', 'lCaud', 'lPut',  'lGP', 'lHipp', 'lAmgd', 'lAcumb', 'lVentDienceph',\
                     'lThirdVent', 'lFourthVent', 'lBrainStem'];
                       
      
      print "::: Converting all ascii segmentations to matlab tri-vert :::";
      for i in range(len(subcort_list)):
          subcort  = subjAscii_dir + '/' + subcort_list[i];
          nuc = subjAscii_dir + '/' + nuc_list[i];
          self.subcortFs2mlab(subcort,nuc);
 

   # method for obtaining freesurfer mni coords using mri_cvs_normalize
   def get_cvsWarp(self):

      '''Method for obtaining freesurfer mni coords using mri_cvs_normalize'''


      # define patient root and fs_dir
      subjects_dir = self.subj_dir;
      patient_root = subjects_dir + self.subj;
      fs_dir = patient_root + '/fs_data/';

      #run cvs register
      orig = fs_dir + 'mri/brain.mgz'; # orig is mri in fs orig space
      template = 'cvs_avg35_inMNI152';
      print "::: Computing Non-linear warping from patient native T1 to fs CVS MNI152 :::";
      os.system('mri_cvs_register --mov %s --template %s &' %(orig,template));
   
      # apply non-linear warp computed from mri_cvs_register
      # to vox coords list to get warped vox and surface coords in MNI
      elecs_dir = '/data_store2/imaging/subjects' + self.subj + '/elecs/';
      cvs_dir = fs_dir + 'cvs/';
      template_brain = subjects_dir + template + '/mri/brain.mgz';
      transform = cvs_dir + 'combined_tocvs_avg35_inMNI152_elreg_afteraseg-norm.tm3d';
      mni_dir = elecs_dir + 'cvs_mni/';
      os.mkdir(mni_dir);
   
      #loop through all patients electrode files to warp each into MNI152
      for coord_list in os.listdir(elecs_dir):
   
         vox_coords = elecs_dir + coord_list;
         mni_voxCoords = elecs_dir + 'mni_' + coord_list;
         print "::: Applying Non-linear warping to %s to obtain %s :::" %(vox_coords, mni_voxCoords);
         os.system('applyMorph --template %s --transform %s point_list %s %s linear' %(template_brain,transform,vox_coords,mni_voxCoords));
      
         # save each as .mat for visualization
         coords = np.loadtxt(mni_voxCoords);
         coords = np.array(coords);
         new_fname = mni_voxCoords.replace('.txt', '.mat');
         scipy.io.savemat(new_fname, {'elecmatrix': coords});

         
                           
   # method calls annot2dpv and splitsrf to tessellate and render meshes 
   # for freesurfer cortical parcellation meshes      
   def get_cortSegMeshes(self):
    
      '''Method calls annot2dpv and splitsrf to tessellate and render meshes
         for freesurfer cortical parcellation meshes'''


      # location of freesurfer look up table for naming ROIs
      fsLUT_dir= self.fsLUT_dir;
    
      # open the LUT text file and grab list of label #'s and names
      fname = fsLUT_dir + '/FreeSurferColorLUT.txt';
      file = open(fname, 'r');
      LUT = file.readlines();
      LUT = LUT[429:502]; # grab only crtx ROIs from .txt LUT
      LUT = [row.split("    ") for row in LUT];
      LUT = [row[0:2] for row in LUT]; # only retain label # and name
      LUT.sort();
      LUT.pop(0); # get rid of non-sense string 
      LUT.pop(0); # and get rid of unknown roi
      LUT.pop(35); # rh unknown
      lh_LUT = LUT[0:35];
      rh_LUT = LUT[35:70];
    
      # nav to subject root dir
      bash_dir = self.bash_dir;
      root = self.subj_dir + '/' + self.subj;
      hem = self.hem;
    
      # convert annot.aparc to .dpv format
      fs_dir = root + '/fs_data';
      os.system('annot2dpv %s/label/%s.aparc.annot %s/label/%s.aparc.annot.dpv' %(fs_dir, hem,fs_dir, hem));
    
      # split and render mesh surface for each cortical roi
      os.system('splitsrf %s/surf/%s.pial.asc %s/label/%s.aparc.annot.dpv %s/surf/%s_pial_roi' %(fs_dir, hem, fs_dir, hem, fs_dir, hem));
    
      # make dir for all pial meshes
      surf_dir = fs_dir + '/surf/';
      os.mkdir('%s/%s_pial_rois' %(surf_dir,hem));
      os.system('mv %s/%s_pial_roi.* %s/%s_pial_rois/' %(surf_dir, hem,surf_dir,hem));
    
      # get list of all the .srf roi files created from splitsrf
      srf_list = [i for i in os.listdir('%s/%s_pial_rois' %(surf_dir,hem)) if (i.endswith('.srf') and i.startswith(hem))];
    
      if hem == 'lh':
        
         # replace roi number with roi name from LUT
         [os.system('mv %s %s' %(i,(lh_LUT[srf_list.index(i)][1] + '.asc'))) for i in srf_list];
       
      elif hem == 'rh':
        
         # replace roi number with roi name from LUT
         [os.system('mv %s %s' %(i,(rh_LUT[srf_list.index(i)][1] + '.asc'))) for i in srf_list]; 
    
      # create corresponding list of ascii file names
      ascii_list = [i for i in os.listdir('%s/%s_pial_rois' %(surf_dir,hem)) if (i.endswith('.asc') and i.startswith('ctx-'+ hem))];
    
      # convert ascii to matlab style tri and vert for meshes
      for fname in ascii_list:
         fname = '%s/%s_pial_rois' %(surf_dir,hem) + '/' + fname;
         self.pialROI_2mlab(fname);
       
      # mkdirs for tri vert and inds
      os.mkdir('%s/%s_pial_rois/tri' %(surf_dir,hem));
      os.mkdir('%s/%s_pial_rois/vert' %(surf_dir,hem));
      os.mkdir('%s/%s_pial_rois/inds' %(surf_dir,hem));
      os.system('mv %s/%s_pial_rois/*tri.mat %s/%s_pial_rois/tri' %(surf_dir,hem,surf_dir,hem));
      os.system('mv %s/%s_pial_rois/*vert.mat %s/%s_pial_rois/vert' %(surf_dir,hem,surf_dir,hem));
      os.system('mv %s/%s_pial_rois/*inds.mat %s/%s_pial_rois/inds' %(surf_dir,hem,surf_dir,hem));
   
      
   # member function for anatomical labeling of electrodes based on euclidean
   # distance from gyrul mesh segmentation (freesurfer)   
   def get_elecLabels(self, elecs):
    
      '''Member function for anatomical labeling of electrodes based on euclidean
         distance from gyrul mesh segmentation (freesurfer). elecs arg is
         a full-path to a .mat file with 'elecmatrix' variable as nx3 
         xyz coordinate list'''

      # 1.) load in elecs and mesh data  
      # nav to subject root dir
      root = self.home_data_dir + self.subj;
      os.chdir(root);
    
      # load electrode coordinates
      os.chdir(root + '/elecs');
      elec_name = elecs.strip('.mat');
      elecs = scipy.io.loadmat(elecs);
      elecs = elecs.get('elecmatrix');
    
      # load pial mesh
      rois_dir = (root + '/fs_data/surf/%s_pial_rois/vert' %(self.hem));
      os.chdir(rois_dir);

      # 2.) form roi list and data containers   
      #create a list of rois to iterate over
      roi_list = [i for i in os.listdir(os.getcwd()) if (i.endswith('.mat') and i.startswith('ctx-'+ self.hem))];
    
      # empty dict for all labeled elecs per roi
      elecs_all = {}; 
    
      # empty dict for all euclidean distances to closest mesh point
      # for each elec (values) per roi (key)
      min_eucDists_all = {};
    
      # dict to store rois as keys, and a list of all electrode indices 
      # for electrodes close in euclidean space to that roi mesh
      roi_elec_inds_all = {};
    
      # convert coords to list. for grabbing indices of electrodes using .index
      elecs_list = list(elecs); # list with each coordinate as an XYZ array
      elecs_list = [list(i) for i in elecs_list];# list with each coordinate as an XYZ list
    
      # iterable list to iterate over each electrode in list
      elec_nums = np.linspace(0,len(elecs)-1, num=(len(elecs)));
 
      # 3.) iterate over roi list to find closest gyrul mesh of each electrode  
      for roi in roi_list:
        
         # get vertex data for current roi
         verts = scipy.io.loadmat(roi);
         verts = verts.get('vert');
    
         # calculate euclidean distance between all electrodes and all mesh points
         eucDist = cdist(elecs, verts, 'euclidean');
    
         # get the minimum euc mesh distance for each electrode and current roi
         min_eucDists = [min(eucDist[i,]) for i in elec_nums];
         min_eucDists = scipy.array(min_eucDists);
    
         # collect all electrodes with a euc distance < 5 to current roi mesh
         roi_elecs = [elecs[i,:] for i in elec_nums if min_eucDists[i] <= (0.55 *scipy.std(min_eucDists))];
         roi_elecs = np.array(list(roi_elecs));
    
         # add selected electrodes to dict entry for current roi
         # elecs_all is a dict with rois as keys and lists of elec coords with 
         # close euclidean distances as values
         elecs_all[roi] = roi_elecs;
       
         # get indices of all electrodes selected for current roi       
         # empty list for poulating with electrode inds
         roi_elec_list = list(roi_elecs);
         roi_elec_list = [list(i) for i in roi_elec_list]; # convert list of arrays to list of lists[xys]
         roi_elec_inds = []; 
         elec_ind_nums = np.linspace(0,len(roi_elecs)-1, num=(len(roi_elecs)));
         elec_ind_nums = elec_ind_nums.astype('int'); # index into roi_elec_list
       
         # populate elec_ind_nums with indeces of all ellectrodes for current roi
         for i in elec_ind_nums:
             roi_elec_inds.append(elecs_list.index(roi_elec_list[i])); # index into original elecs_lsit
       
         # add roi elec inds to dict of all roi elec inds: each roi(key) has a list
         # of indices(val)
         roi_elec_inds_all[roi] = roi_elec_inds;
       
         # append euclidian distances to current mesh for all electrodes
         # min_euc_Dists_all is a dict with rois for keys, each with a list of
         # euclidean distances from each electrode (for the closest point on the
         # current mesh) 
         min_eucDists_all[roi] = min_eucDists;
       

      # 4.) find electrodes that are close to two rois and compare to find closest 
      # roi for a given electrode 
  
      # closest euclidean distance for an electrode
      # get dict of electrode indices (w), with values being a list of rois overlapping
      w = collections.defaultdict(list); 
      for k,v, in roi_elec_inds_all.iteritems():
         for i in v: w[i].append(k);
    
      # iterate over electrode indexes
      for k in w.keys():
           
         # if for a given electrode index, there are two or more rois listed with the same electrode
         if len(w[k]) >1:
            min_dist = 1000000; # set initial min dist to 1 mil
             
            # find which region has a smaller euc distance to that electrode
            for region in w[k]:
                 
               # if the current region has smaller eucdist make it the min
               if min_eucDists_all[region][k] < min_dist:
                  min_dist = min_eucDists_all[region][k];
                  closest = region;
                     
               # change the key label to the closest region
               w[k] = closest;
                 

      # output final list of each electrodes label
      elec_labels =[];
      for k in w.keys():
          if type(w[k]) == list:
             elec_labels.append(w[k][0]);
          else:
             elec_labels.append(w[k]);
            
      # get rid of 'ctx-' and '_vert.mat' for each label
      elec_labels = [label.strip('_vert.mat') for label in elec_labels];
      elec_labels = [label.strip('ctx-') for label in elec_labels];
    
      # list of all rois included in final labeling
      label_rois = set(elec_labels);
      label_rois = list(label_rois);
        
    
      # 5.) Save final outputs  
      # create dirs for saved labeled elecs
      os.chdir(root + '/elecs');
      elecs_dir = os.getcwd();
      labels_dir = elecs_dir + '/labels';
    
      # check if labels dir already exists
      if os.path.isdir(labels_dir) == True:
        
          None;
        
      else:
        
          os.mkdir('labels');
    
      # go in to labels dir   
      os.chdir(labels_dir);
    
      # save coords .mat for each roi
      elec_nums = [int(num) for num in elec_nums];
    
      #check if any electrodes went missing and didn't get labeled
      elec_nums_set = set(elec_nums);
      missing = elec_nums_set.difference(w.keys());
      missing = list(missing);
    
      # if there are missing electrodes add the missing electrode indeces to
      # the list of electrode keys, and fill in the label using the adjacent 
      # electrode labels in elec_labels
      if missing != None:
         keys = w.keys()
         for missed in missing:
             keys.insert(missed, missed);
             elec_labels.insert(missed,elec_labels[missed]);
           
    
         keys.sort(); # make sure keys is in order
          
           
        
      # save each list of labeld elecs
      for roi in label_rois:
        
          # container for elecs on current roi to save
          labeled_elecs = [];
        
          # add all electrodes into container if label is current roi
          for i in elec_nums:
              if elec_labels[i] == roi:
                  labeled_elecs.append(elecs[i]);
        
          # save a .mat for the current roi coords
          labeled_elecs = np.array(labeled_elecs);
          fname = elec_name + '_' + roi + '_elecs.mat';
          scipy.io.savemat(fname, {'elecmatrix':labeled_elecs});


      # return list of elec_labels
      return elec_labels;
      
      
      
   # member function for cleaning outputs of 3d spherical
   # hough electrode detection.    
   def clean_Hough3D(self,fname):

      '''Member function for cleaning outputs of 3d spherical
         hough electrode detection. fname arg is a full-path
         to a .mat file containing a variable 'sphcen2' w/
         xyz matrix (nx3) of elec center coordinates'''


      # load in 3D hough output
      coords = scipy.io.loadmat(self.CT_dir + '/' + fname).get('sphcen2');
      x = coords[:,0];
      y = coords[:,1];
      z = coords[:,2];

      # plot dirty elec results
      s = 150;
      fig = plt.figure(facecolor="white");
      ax = fig.add_subplot(111,projection='3d');
    

      # create point and click function for selecting false alarm points
      false_elecs = [];
      def onpick(event):

         # get event indices and x y z coord for click
         ind = event.ind[0];
         x_e, y_e, z_e = event.artist._offsets3d;
         print x_e[ind], y_e[ind], z_e[ind];
         false_coord = [x_e[ind], y_e[ind], z_e[ind]];
         false_elecs.append(false_coord);

         # refresh plot with red dots picked and blue dots clean
         #fig.clf();
         false_arr = np.array(false_elecs);
         x_f = false_arr[:,0];
         y_f = false_arr[:,1];
         z_f = false_arr[:,2];
         plt.cla();
         Axes3D.scatter3D(ax,x_f,y_f,z_f,s=150,c='r', picker=5);

         # get array of only clean elecs to re-plot as blue dots
         clean_list = list(coords);
         clean_list = [list(i) for i in clean_list];
         for coordin in clean_list:
                if list(coordin) in false_elecs:
                    clean_list.pop(clean_list.index(coordin));
         clean_coordis = np.array(clean_list);
         x_c = clean_coordis[:,0];
         y_c = clean_coordis[:,1];
         z_c = clean_coordis[:,2];
         Axes3D.scatter3D(ax,x_c, y_c, z_c, s=150,c='b', picker=5);
         time.sleep(0.5);
         plt.draw();
    
      Axes3D.scatter3D(ax,x,y,z,s=s, picker=5);
      #title_font = {'fontname':'serif', 'size':'24', 'color':'black', 'weight':'bold',                            
      #'verticalalignment':'bottom'};                                                                        
      plt.title('Hough 3D Output Coordinates: ');
      fig.canvas.mpl_connect('pick_event', onpick);
      plt.ion();
      plt.show(fig);

      # wait for user to press enter to continue
      while True:
        
         fig.canvas.mpl_connect('pick_event', onpick);
         plt.draw();
         i = raw_input("Enter text or Enter to quit: ");
         if i == '':
            
             # remove all false coords from elec centers list
             fig.clf();
             coords_list = list(coords);
             coords_list = [list(i) for i in coords_list];
             for coordi in coords_list:
                 if list(coordi) in false_elecs:
                     coords_list.pop(coords_list.index(coordi));
             clean_coords = np.array(coords_list)

             # plot result
             X = clean_coords[:,0];
             Y = clean_coords[:,1];
             Z = clean_coords[:,2];


             s = 150;
             fig = plt.figure(facecolor="white");
             ax = fig.add_subplot(111,projection='3d');
             Axes3D.scatter3D(ax,X,Y,Z,s=s);
             plt.show(fig);
            
             # save the cleaned elecs file
             new_fname = self.CT_dir + '/' + 'Hough_3D_cleaned.mat';            
             scipy.io.savemat(new_fname, {'elecmatrix':clean_coords});
            
             return clean_coords;
             break;



   # member function for windowed intensity based thresholding of an img
   # thresh_f is the intensity floor, and thresh_c is the intensity ceiling
   def electhresh(self,img,thresh_f,thresh_c):

      '''Member function for windowed intensity based thresholding of an img
         thresh_f is the intensity floor, and thresh_c is the intensity ceiling'''
    
      # function thresholds a CT to pinpoint electrode locations
      img_fname = img;
      img = nib.load(img);
      img_dat = img.get_data();
      img_arr = np.array(img_dat);
      new_fname = 'thresh_' + img_fname;
      img.set_filename = new_fname;
    
      # Create a container for the thresholded CT values
      CT_thresh = np.zeros(img_arr.shape);
    
      # Get number of slices in CT
      Z_dim = np.size(img_arr[1,1,:]);
      X_dim = np.size(img_arr[:,1,1]);
      y_dim = np.size(img_arr[1,:,1]);
    
      # Threshold the CT
      for z in range(1,Z_dim):
         CT_thresh[:,:,z] = ((img_arr[:,:,z] > (thresh_f * img_arr[:,:, z].std())) & (img_arr[:,:,z] < (thresh_c * img_arr[:,:, z].std())));
         for x in range(1,X_dim):
            for y in range(1,y_dim):
               if CT_thresh[x,y,z] > 0:
                  CT_thresh[x,y,z] = 1;
                
      # save the thresholded ct
      hdr = img.get_header();
      affine = img.get_affine();
      N = nib.Nifti1Image(CT_thresh,affine,hdr);
      N.to_filename(new_fname);
      

   # member function for interpolating depth probe
   # contact coordinates from two endpoint coordinates
   def interp_depth(self,elec1, elec2, numElecs, fname):

      '''Member function for interpolating depth probe
         contact coordinates from two endpoint coordinates.
         elec1 is a full-path to a .mat file with xyz
         coord of first electrode of one end of a depth probe.
         elec2 is the last end coordinate. numElecs arg is
         the number of electrodes on the probe total. 
         fname is the name you want to save the output to.'''

      # load in elec1 and elec2 coords
      elec1 = scipy.io.loadmat(elec1).get('elecmatrix');
      elec2 = scipy.io.loadmat(elec1).get('elecmatrix');

      # interpolate through x y and z from end to end
      xs = np.linspace(elec1[0], elec2[0], numElecs);
      ys = np.linspace(elec1[1], elec2[1], numElecs);
      zs = np.linspace(elec1[2], elec2[2], numElecs);

      # flip interp xyz coords to stack as one matrix
      xs = np.transpose(xs);
      ys = np.transpose(ys);
      zs = np.transpose(zs);
      stackem = (xs,ys,zs);
      elecmatrix = np.column_stack(stackem);

      #  save result
      scipy.io.savemat(fname, {'elecmatrix':elecmatrix});
      
      
      
   # member function for converting a .mat electrode coordinate list
   # to image coordinates in a NIFTI volume for quality checking
   def elecs2brainVox(self, img, elecs_file, voxCoords_fname, elecsVol_fname):

      '''Member function for converting a .mat electrode coordinate list
         to image coordinates in a NIFTI volume for quality checking'''

      # load fs brain.nii and get image data
      nat_brain = nib.load(img);
      dat = nat_brain.get_data();

      # create empty image to put electrodes in 
      elecs_img = numpy.zeros(dat.shape);

      # chdir to elec coordinates file
      elecs = scipy.io.loadmat(elecs_file)
      elecs = elecs.get('elecmatrix');

      # vox2Ras affine for all fs orig volumes
      affine = numpy.array([[  -1.,    0.,    0.,  128.],
             [   0.,    0.,    1., -128.],
             [   0.,   -1.,    0.,  128.],
             [   0.,    0.,    0.,    1.]]);
       
      # get inverse of orig vox2ras
      inv_affine = numpy.linalg.inv(affine);

      # apply invers RAS2vox to electrode coords
      intercept = numpy.ones(len(elecs));      
      vox_elecs = numpy.column_stack((elecs,intercept));
      vox_elecs = numpy.transpose(vox_elecs);
      elec_vox = numpy.matrix.dot(inv_affine,vox_elecs);
      elec_vox = numpy.transpose(elec_vox);
      elec_vox = elec_vox[:,0:3];
      elec_vox = elec_vox.astype('int');

      # put ones in the vox coords of elecs in empty img
      for i in elec_vox:
         elecs_img[i[0], i[1], i[2]] = 1;

      # write vox image to nifti
      hdr = nat_brain.get_header();
      affine = nat_brain.get_affine();
      fname = elecsVol_fname;
      N = nib.Nifti1Image(elecs_img,affine,hdr);
      N.to_filename(fname);
   
      # save as .txt
      numpy.savetxt(voxCoords_fname, elec_vox, delimiter=' ', fmt='%-7.1f');
   
      # return electrode voxel coordinates in fs orig space
      return elec_vox;
      
      
   # member function to resample a high density grid (256) 
   # to clinical low-desnity montage (64)  
   def get_clinGrid(self):
    
       '''member function to resample a high density grid (256)
          to clinical low-desnity montage (64)'''

       # define elecs dir
       elecs_dir = self.elecs_dir;
    
       # nav to elecs with elecs file
       os.chdir(elecs_dir);
    
       # inds of electrodes recorded in clinical montage
       inds = np.array([[  0,   2,   4,   6,   8,  10,  12,  14,  32,  34,  36,  38,  40,
            42,  44,  46,  64,  66,  68,  70,  72,  74,  76,  78,  96,  98,
           100, 102, 104, 106, 108, 110, 128, 130, 132, 134, 136, 138, 140,
           142, 160, 162, 164, 166, 168, 170, 172, 174, 192, 194, 196, 198,
           200, 202, 204, 206, 224, 226, 228, 230, 232, 234, 236, 238]]);

       # load high density grid
       hd = scipy.io.loadmat('hd_grid.mat');
       hd = hd.get('elecmatrix');

       # resample using indices of electrode coordinates
       clinGrid = hd[inds,:];
       clinGrid = clinGrid.reshape(64,3); #get rid of erroneuos third dimension

       # save resmpled clinical grid
       scipy.io.savemat('clinical_grid.mat', {'elecmatrix':clinGrid});  
       
       
       
   # function for running streamlines tractography between 
   # all freesurfer rois for a given subject. 
   def get_Connectome(self):
    

      # set up var for DWI dir and roi dir
      DWI_dir = self.DWI_dir;
      roi_dir = self.subj_dir + self.subj + '/mri/crtx_subcort_ROIs/';
      os.chdir(roi_dir);
   
      # get list of all rois
      roi_list = os.listdir(os.getcwd());
      roi_list2 = os.listdir(os.getcwd());
   
      # create dir for ouput of tractography
      os.mkdir('roi_roi_tcks');

   
      # run roi-roi tractography for all roi pars
      for roi in roi_list:
          for roi_s in roi_list2:

              roi1 = roi_dir + '/' + roi;
              roi2 = roi_dir + '/' + roi_s;
              command = 'tckgen %sFOD.mif %s_%s.tck -seed_image %s -include %s -mask \
                      %smask.mif -nthreads 12 -maxnum 10000' %(self.DWI_dir, roi, roi_s, roi1, roi2, self.DWI_dir);
              args = shlex.split(command);
              subprocess.Popen(args);
    
      exit# clean up results              
      os.system('mv *.tck roi_roi_tcks/');
      
      
      
   # function gets freesurfer c_RAS vector
   # to save for shifting meshes into orig space correctly
   def get_cRAS(self):
    
      '''Function gets freesurfer c_RAS vector 
         to save for shifting meshes into orig space correctly.'''

      # define subject dir
      subj_fs_dir = self.subj_dir;
    
      # use mris_info to get c_RAS from lh.pial freesurfer ouput
      os.system('mris_info %s/surf/lh.pial >& %s/surf/lh_pial_info.txt' %(subj_fs_dir, subj_fs_dir));
      command = 'grep "c_(ras)" %s/surf/lh_pial_info.txt' %(subj_fs_dir);
      args = shlex.split(command);
      p = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0];
      p = p.strip('\n'); p = p.strip('c_(ras): '); p = p.strip('(');p = p.strip(')');
      p = p.split(','); p = [float(i) for i in p]; p = np.array(p);
    
      c_RAS = p; # this is the coordinate shift from lh.pial to lh.pial.asc coords
      scipy.io.savemat('%s/surf/c_RAS.mat' %(self.subj_dir), {'c_RAS':c_RAS});
    
      return c_RAS; 


   # function to clean outputs of AFQ wm bundle segmentation on CSD data
   def clean_fiberOutliers(self):
       
       mlab = matlab.MatlabCommand();
       mlab.inputs.script = "addpath(genpath('%s')); \
                             addpath(genpath('%s')); \
                             addpath(genpath('%s')); \
                             load('%s/dti55trilin/CSD_segmented_fibers.mat'); \
                             maxDist = 2.6; maxLen = 100; \
                             numNodes = 25; M = 'mean'; count = 1; show = 0; \
                             for i = 12:20 \
                                fg = fgs(i); \
                                [fgclean keep]=AFQ_removeFiberOutliers(fg,maxDist,maxLen,numNodes,M,count,show); \
                                fgs(i) = fgclean; \
                             end; \
                             save('%s/dti55trilin/CSD_segmented_fibers_clean.mat', 'fgs');" %(self.AFQ_dir, self.SPM_dir, self.VISTA_dir, self.DWI_dir, self.DWI_dir);
       
       
       print "::: Cleaning AFQ WM Bundle Segmentation :::";        
       out = mlab.run();
       print "::: Finished :::"

       return out;
          
          
   # script takes dtiInit output DWI and produces masked DWI and binary mask
   def get_DWImask(self):
    
      #nav to DWI dir
      os.chdir(self.DWI_dir);
    
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




   # Script extracts tensor models, FA, and major eigenvectors from HARDI diffusion 
   # data pre-processed with vistasoft-AFQ dtiInit pipeline (disotrtion-eddy current 
   # correction), and then extracts a response function from the DWI, and uses it to 
   # perform constrained spherical deconvolution and probobailistic fiber tracking
   # on the resultant FOD fields
   def get_mrtrix_metrics(self):
    
      # nav to DWI dir
      os.chdir(self.DWI_dir);
    
      # use dipy to pull out white matter mask
      self.get_DWImask(self.DWI_dir);
    
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
      
      
   # member function to apply a guassian smoth to a NIFTI img
   def img_filt(self, fname):
    
      # Load the thresholded CT
      img = nib.load(fname);

      # convert to nipy format for filtering
      nipy_img = nipy.io.nifti_ref.nifti2nipy(img);

      # build and apply gaussian filter
      smoother = LinearFilter(nipy_img.coordmap, nipy_img.shape, 1);
      smoothed_im = smoother.smooth(nipy_img);

      # convert smoothed im back to nifti
      nifti_img = nipy.io.nifti_ref.nipy2nifti(smoothed_im);

      # get data array for saving
      smoothed_CT = nifti_img.get_data();
      smoothed_CT = np.array(smoothed_CT);

      # save the smoothed thresholded ct
      hdr = img.get_header();
      affine = img.get_affine();
      N = nib.Nifti1Image(smoothed_CT,affine,hdr);
      new_fname = 'smoothed_' + fname;
      N.to_filename(new_fname);
 
    
   def make_CortROIs(self):
    
       # location of freesurfer look up table for naming ROIs
       fsLUT_dir= self.fsLUT_dir;
       os.chdir(fsLUT_dir);
    
       # open the LUT text file and grab list of label #'s and names
       fname = 'FreeSurferColorLUT.txt';
       file = open(fname, 'r');
       LUT = file.readlines();
       cort_LUT = LUT[429:502]; # grab only crtx ROIs from .txt LUT
       cort_LUT = [row.split("    ") for row in cort_LUT];
       cort_LUT = [row[0:2] for row in cort_LUT]; # only retain label # and name
       cort_LUT.sort();
       cort_LUT.pop(0); # get rid of non-sense string
    
       # get subcort seg lists too
       lh_subcort = [['8', ' Left-Cerebellum-Cortex'],
       ['9', ' Left-Thalamus'],
       ['10', 'Left-Thalamus-Proper'],
       ['11', 'Left-Caudate'],
       ['12', 'Left-Putamen'],
       ['13', 'Left-Pallidum'],
       ['16', 'Brain-Stem'],
       ['17', 'Left-Hippocampus'],
       ['18', 'Left-Amygdala'],
       ['19', 'Left-Insula'],
       ['20', 'Left-Operculum'],
       ['26', 'Left-Accumbens-area'],
       ['27', 'Left-Substancia-Nigra'],
       ['28', 'Left-VentralDC']];
    
       # r hem subcort
       rh_subcort = [['47', 'Right-Cerebellum-Cortex'],
       ['48', 'Right-Thalamus'],
       ['49', 'Right-Thalamus-Proper'],
       ['50', 'Right-Caudate'],
       ['51', 'Right-Putamen'],
       ['52', 'Right-Pallidum'],
       ['53', 'Right-Hippocampus'],
       ['54', 'Right-Amygdala'],
       ['55', 'Right-Insula'],
       ['56', 'Right-Operculum'],
       ['58', 'Right-Accumbens-area'],
       ['59', 'Right-Substancia-Nigra'],
       ['60', 'Right-VentralDC']];
    
       # once list with all roi #s and names
       roi_list = cort_LUT + lh_subcort + rh_subcort;

       # nav to patient's fs_data/mri dir 
       imgs_dir = '/Volumes/Thunderbolt_Duo/imaging/subjects/';
       patient_dir = imgs_dir + self.subj + '/fs_data/mri';
       os.chdir(patient_dir);
    
       # convert aparc+aseg to NIFTI
       os.system('mri_convert aparc+aseg.mgz aparc_aseg.nii');
    
       # loop through LUT to make roi for each label
       for row in roi_list:
           label_num = int(row[0]);
           label_name = row[1];
           print "::: Creating ROI for %s :::" %(label_name);
           command1 = 'fslmaths aparc_aseg.nii -uthr %d -thr %d %s.nii' %(label_num, label_num, label_name);# create label ROI
           args1 = shlex.split(command1);
           p1 = subprocess.Popen(args1);
           command2 = 'fslmaths %s.nii.gz -div %d %s.nii.gz' %(label_name, label_num, label_name); # binarize mask to 1s and 0s
           args2 = shlex.split(command2);
           p2 = subprocess.Popen(args1);
        
       # make dir to move all labels into
       os.mkdir('crtx_subcort_ROIs');
       os.system('mv *.nii.gz crtx_subcort_ROIs/');
    
       print "::: Finished :::";   

   # member function to plot brain surface as triangular mesh
   def ctmr_brain_plot(self,tri,vert):
    
      # load in mesh data
      vert = scipy.io.loadmat(vert).get('vert')
      tri = scipy.io.loadmat(tri).get('tri');
      tri = tri -1; # correct tri from matlab style index
    
      # define x y z coords of vertices
      x = vert[:,0];
      y = vert[:,1];
      z = vert[:,2];
    
      # plot cortex and begin display
      mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1));
      mesh = mlab.triangular_mesh(x,y,z,tri, color=(0.8,0.8,0.8), resolution=200);

      # change OpenGL mesh properties for phong point light shading
      mesh.actor.property.ambient = 0.0;
      mesh.actor.property.specular = 0.2353;
      mesh.actor.property.specular_power = 128;
      mesh.actor.property.diffuse = 1.0;
      mesh.actor.property.interpolation = 'phong';
      mesh.scene.light_manager.light_mode = 'vtk';


      return mesh;


   # member function to plot electrode coordinates
   def el_add(self, elecs, color, s):

      #load in elecs coords
      elecs = scipy.io.loadmat(elecs).get('elecmatrix');

      # define x y z for mayavi sphere plot
      x_elec = elecs[:,0];
      y_elec = elecs[:,1];
      z_elec = elecs[:,2];

      # plot the electrodes as spheres
      points = mlab.points3d(x_elec,y_elec,z_elec, scale_factor = s, color = color, resolution=25);
      points.actor.property.ambient = 0.049;
      points.actor.property.specular = 1;
      points.actor.property.specular_power = 66;
      points.actor.property.diffuse = 0.7059;
      points.actor.property.interpolation = 'phong';
      points.scene.light_manager.light_mode = 'vtk';

      return points;

   # plot guasian smoothed colored surf using scalar data
   def paint_brain(self, tri,vert,elecs,wts, guass):

       # load in mesh data                                                                                       
       vert = scipy.io.loadmat(vert).get('vert')
       tri = scipy.io.loadmat(tri).get('tri');
       tri = tri -1; # correct tri from matlab style index                                                       

       # define x y z coords of vertices                                                                         
       x = vert[:,0];
       y = vert[:,1];
       z = vert[:,2];

       # load in elecs
       elecs = scipy.io.loadmat(elecs).get('elecmatrix');

       # find closest vertices to electrodes and set scalar value at those verts
       # according to the input wts
       s = np.zeros(len(x))
       eucDist = cdist(elecs, vert, 'euclidean');
       elec_nums = np.linspace(0,len(elecs)-1, num=(len(elecs)));
       elec_vert = [np.where(row == row.min()) for row in eucDist]
       elec_vert = scipy.array(elec_vert)
       elec_vert = elec_vert[:,0,0]

       # assign the scalar of the found closest vertices to appropriate wts
       s_list = list(s);
       elec_vert_list = list(elec_vert);
       for v in elec_vert_list:
          index = elec_vert_list.index(v);
          s[v] = wts[index] * guass;


       # plot cortex and begin display                                                                           
       mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1));
       mesh = mlab.triangular_mesh(x,y,z,tri, scalars=s, colormap='RdBu');

       # change OpenGL mesh properties for phong point light shading                                             
       mesh.actor.property.ambient = 0.0;
       mesh.actor.property.specular = 0.2353;
       mesh.actor.property.specular_power = 128;
       mesh.actor.property.diffuse = 1.0;
       mesh.actor.property.interpolation = 'phong';
       mesh.scene.light_manager.light_mode = 'vtk';

       return mesh;

   # plot AFQ style DWI tractography bundle as tube
   def plot_bundle(self, Bundle,Color,tube_size):

      # load in DWI bundle 
      bundle = scipy.io.loadmat(Bundle).get('fg');
      num_fibers = len(bundle[0,0][8]);

      # plot the DWI streamline fibers as tubes
      for fib in range(num_fibers):
          fiber = bundle[0,0][8][fib][0];
          fiber = np.transpose(fiber);
          fiber = np.array(fiber);
          x = fiber[:,0];
          y = fiber[:,1];
          z = fiber[:,2];
          fiber_mesh = mlab.plot3d(x,y,z, color=Color, tube_sides = 33, tube_radius=tube_size);

          # change mesh properties for shiny tubes
          fiber_mesh.actor.property.ambient = 0.0;
          fiber_mesh.actor.property.specular = 0.2353;
          fiber_mesh.actor.property.specular_power = 128;
          fiber_mesh.actor.property.diffuse = 1.0;
          fiber_mesh.actor.property.interpolation = 'phong';
          fiber_mesh.scene.light_manager.light_mode = 'vtk';
          
          
   # plot translucent glass brain surface 
   def glass_brain(self, tri,vert):
    
      # load in mesh data
      vert = scipy.io.loadmat(vert).get('vert')
      tri = scipy.io.loadmat(tri).get('tri');
      tri = tri -1; # correct tri from matlab style index
    
      # define x y z coords of vertices
      x = vert[:,0];
      y = vert[:,1];
      z = vert[:,2];
    
      # plot cortex and begin display
      mlab.figure(fgcolor=(0, 0, 0), bgcolor=(0, 0, 0));
      mesh = mlab.triangular_mesh(x,y,z,tri, color=(0.8,0.8,0.8), resolution=200);

      # change OpenGL mesh properties for phong point light shading
      mesh.actor.property.ambient = 0.186;
      mesh.actor.property.specular = 1;
      mesh.actor.property.specular_power = 128;
      mesh.actor.property.diffuse = 0.186;
      mesh.actor.property.opacity = 0.2481;
      mesh.actor.property.interpolation = 'phong';
      mesh.scene.light_manager.light_mode = 'vtk';


      return mesh;
      

   # function to run DWI pre-processing and way-point roi
   # based tractography of major white-matter bundles
   def DWI_preproc(self):
       
       mlab = matlab.MatlabCommand();
       mlab.inputs.script = "addpath(genpath('%s')); \
                             addpath(genpath('%s')); \
                             addpath(genpath('%s')); \
                             raw_DWI_dir = '%s'; \
                             raw_DWI_dir = [raw_DWI_dir '/raw/']; \
                             dwParams = dtiInitParams; \
                             dwParams.bvecsFile = [raw_DWI_dir 'dwi.bvec']; \
                             dwParams.bvalsFile = [raw_DWI_dir 'dwi.bval']; \
                             dwParams.phaseEncodeDir = 2; \
                             dwParams.bvalue = %d; \
                             dwRaw = [raw_DWI_dir 'dwi.nii.gz']; \
                             [dt6FileName, outBaseDir] = dtiInit(dwRaw,[raw_DWI_dir 'orig_brain.nii'],dwParams); \
                             dt = dtiLoadDt6(['%s/dti55trilin/dt6.mat']); \
                             fg = AFQ_WholebrainTractography(dt); \
                             save([raw_DWI_dir 'dti55trilin/DTI_fg_wholeBrain_tracts.mat'], 'fg'); \
                             fgs = AFQ_SegmentFiberGroups(dt); \
                             save([raw_DWI_Dir 'dti55trilin/DTI_segmented_fibers.mat'], 'fgs');" %(self.AFQ_dir, self.SPM_dir, self.VISTA_dir, self.DWI_dir, self.DWI_bval, self.DWI_dir);
    
       print "::: Initiating DWI pre-processing and WM Bundle Segmentation :::";        
       out = mlab.run();
       print "::: Finished :::";
       
       
       
   # function to copy fs brain.mgz to DWI dir
   # and convert both DWI and brain to LAS orientation before pre-proc
   def prep_DWI_brain(self):
        
      # convert fs brain.mgz to .nii
      fs_brain = self.subj_dir + self.subj + '/mri/brain.mgz';
      new_fname = fs_brain.replace('brain.mgz', 'orig_brain.nii')
      os.system('mri_convert %s %s' %(fs_brain, new_fname));
      self.fs_brain = fs_brain;
      
      # copy fs brain.nii to DWI raw dir
      DWI_raw_dir = self.DWI_dir + '/raw/';
      self.DWI_raw_dir = DWI_raw_dir;
      os.system('cp %s %s' %(new_fname, DWI_raw_dir));
      
      # reorient both brain.nii and DWI.nii to LAS orientation
      self.DWI_brain = self.DWI_raw_dir + '/orig_brain.nii';
      os.system('fslorient -forceneurological %s' %(self.DWI_brain));
      self.DWI_raw = self.DWI_raw_dir + '/dwi.nii.gz';
      os.system('fslorient -forceneurological %s' %(self.DWI_raw));
      
   
   # function to convert freesurfer ascii subcort segmentations
   # to triangular mesh array .mat style   
   def subcortFs2mlab(self, subcort,nuc):
    
      '''Function to convert freesurfer ascii subcort segmentations
         to triangular mesh array .mat style.'''

      # use freesurfer mris_convert to get ascii subcortical surface
      subcort_ascii = subcort;
    
      # clean up ascii file and extract matrix dimensions from header
      subcort = open(subcort_ascii, 'r');
      subcort_mat = subcort.readlines();
      subcort.close();
      subcort_mat.pop(0); # get rid of comments in header
      subcort_mat = [item.strip('\n') for item in subcort_mat];# get rid of new line char
    
      # extraxt inds for vert and tri    
      subcort_inds = subcort_mat.pop(0); 
      subcort_inds = subcort_inds.split(' '); # seperate inds into two strings
      subcort_inds = [int(i) for i in subcort_inds];# convert string to ints
    

      # get rows for vertices only, strip 0 column, and split into seperate strings
      subcort_vert = [item.strip(' 0') for item in subcort_mat[:subcort_inds[0]]];
      subcort_vert = [item.split('  ') for item in subcort_vert]; # seperate strings
    
      # create containers for each column in vert matrix
      x = [];
      y = [];
      z = [];
    
      # fill the containers with float values of the strings in each column
      for i in subcort_vert:
          x.append(float(i[0]));
          y.append(float(i[1]));
          z.append(float(i[2]));
        
      # convert to scipy mat
      x = scipy.mat(x);
      y = scipy.mat(y);
      z = scipy.mat(z);
    
      # concat columns to one n x 3 matrix
      x = x.transpose();
      y = y.transpose();
      z = z.transpose();
      subcort_vert = scipy.column_stack((x,y,z));
      scipy.io.savemat('subcort_%s_vert.mat' %(nuc), {'vert':subcort_vert});# save vert matrix
    
      # get rows for triangles only, strip 0 column, and split into seperate strings
      subcort_tri = [item[:-2] for item in subcort_mat[subcort_inds[0] +1:]];
      subcort_tri = [item.split(' ') for item in subcort_tri]; # seperate strings  
    
      # create containers for each column in vert matrix
      x = [];
      y = [];
      z = [];
    
      # fill the containers with float values of the strings in each column
      for i in subcort_tri:
          x.append(int(i[0]));
          y.append(int(i[1]));
          z.append(int(i[2]));
        
      # convert to scipy mat
      x = scipy.mat(x);
      y = scipy.mat(y);
      z = scipy.mat(z);
    
      # concat columns to one n x 3 matrix
      x = x.transpose();
      y = y.transpose();
      z = z.transpose();
      subcort_tri = scipy.column_stack((x,y,z));
      subcort_tri = subcort_tri + 1;
      scipy.io.savemat('subcort_%s_tri.mat' %(nuc), {'tri':subcort_tri});# save tri matrix  

      # convert inds to scipy mat
      subcort_inds = scipy.mat(subcort_inds);
      scipy.io.savemat('subcort__%s_inds.mat' %(nuc), {'inds':subcort_inds});# save inds   
      

   # function to convert pial roi segmentations from freesurfer
   # to matlab style .mat triangular mesh data
   def pialROI_2mlab(self, fname):
    
       '''Function to convert pial roi segmentations from freesurfer
          to matlab style .mat triangular mesh data.'''

       # use freesurfer mris_convert to get ascii subcortical surface
       roi_ascii = fname;
    
       # clean up ascii file and extract matrix dimensions from header
       ROI = open(roi_ascii, 'r');
       roi_mat = ROI.readlines();
       ROI.close();
       roi_mat.pop(0); # get rid of comments in header
       roi_mat = [item.strip('\n') for item in roi_mat];# get rid of new line char
    
       # extraxt inds for vert and tri    
       roi_inds = roi_mat.pop(0); 
       roi_inds = roi_inds.split(' '); # seperate inds into two strings
       roi_inds = [int(i) for i in roi_inds];# convert string to ints
      

       # get rows for vertices only, strip 0 column, and split into seperate strings
       roi_vert = [item.strip(' 0') for item in roi_mat[:roi_inds[0]]];
       roi_vert = [item.split(' ') for item in roi_vert]; # seperate strings
    
       # create containers for each column in vert matrix
       x = [];
       y = [];
       z = [];
    
       # fill the containers with float values of the strings in each column
       for i in roi_vert:
           x.append(float(i[0]));
           y.append(float(i[1]));
           z.append(float(i[2]));
        
       # convert to scipy mat
       x = scipy.mat(x);
       y = scipy.mat(y);
       z = scipy.mat(z);
    
       # concat columns to one n x 3 matrix
       x = x.transpose();
       y = y.transpose();
       z = z.transpose();
       roi_vert = scipy.column_stack((x,y,z));
       new_fname = fname.replace('.asc', '_vert.mat');
       scipy.io.savemat(new_fname, {'vert':roi_vert});# save vert matrix
    
       # get rows for triangles only, strip 0 column, and split into seperate strings
       roi_tri = [item[:-2] for item in roi_mat[roi_inds[0] +1:]];
       roi_tri = [item.split(' ') for item in roi_tri]; # seperate strings  
    
       # create containers for each column in vert matrix
       x = [];
       y = [];
       z = [];
    
       # fill the containers with float values of the strings in each column
       for i in roi_tri:
           x.append(int(i[0]));
           y.append(int(i[1]));
           z.append(int(i[2]));
        
       # convert to scipy mat
       x = scipy.mat(x);
       y = scipy.mat(y);
       z = scipy.mat(z);
    
       # concat columns to one n x 3 matrix
       x = x.transpose();
       y = y.transpose();
       z = z.transpose();
       roi_tri = scipy.column_stack((x,y,z));
       roi_tri = roi_tri + 1;
       new_fname = fname.replace('.asc', '_tri.mat');
       scipy.io.savemat(new_fname, {'tri':roi_tri});# save tri matrix  

       # convert inds to scipy mat
       roi_inds = scipy.mat(roi_inds);
       new_fname = fname.replace('.asc', '_inds.mat');
       scipy.io.savemat(new_fname, {'inds':roi_inds});# save inds   
       

   # member function to compute spherical hough transform for localizing electrode centers from post-implant CT
   def Spherical_Hough(self):
       
       '''Member function for computing spherical hough transform for 
          localizing electrode centers from post-implant CT.'''

       # compute hough transform to find electrode center coordinates
       mlab = matlab.MatlabCommand();
       mlab.inputs.script = "addpath(genpath('%s')); addpath(genpath('%s')); addpath(genpath('%s')); addpath(genpath('%s')); \
       img = readFileNifti('%s/smoothed_thresh_CT.nii');img = img.data(); \
       fltr4LM_R = 3; obj_cint = 0.0; radrange = [0.1, 3]; multirad = 0.3; \
       grdthres = 0.33; \
       [center_img,sphere_img,sphcen2,sphrad]=SphericalHough(img,radrange,grdthres, \
       fltr4LM_R, multirad, obj_cint); save('%s/elec_spheres.mat', 'sphere_img'); \
       save('%s/elec_centers.mat', 'center_img'); \
       save('%s/3d_hough_elect_coords.mat', 'sphcen2');" %(self.Hough_dir, self.VISTA_dir, self.SPM_dir, self.VISTA_dir, self.CT_dir, self.CT_dir , self.CT_dir, self.CT_dir);
    
       print "::: Applying 3D hough transform to smoothed thresh-CT :::";        
       out = mlab.run(); # run the hough and get errors if any ocur

       # save txt version of hough putputs
       coords = scipy.io.loadmat(self.CT_dir + '/'+ '3d_hough_elect_coords.mat').get('sphcen2')
       np.savetxt((self.CT_dir + '/dirty_elecsAll.txt'), coords, delimiter=' ', fmt='%-7.1f');
       
       # load native CT to get affine and header
       fname = self.CT_dir + '/smoothed_thresh_CT.nii'
       img = nib.load(fname);
       hdr = img.get_header();
       affine = img.get_affine();

       # save electrode center voxel coordinates to NIFTI img for later use
       img_data = scipy.io.loadmat(self.CT_dir + '/elec_centers.mat')
       img_data = img_data.get('center_img')
       new_fname = self.CT_dir + '/CT_elecCenters.nii'
       N = nib.Nifti1Image(img_data,affine,hdr);
       N.to_filename(new_fname);


       return coords, out;


   # member function to obtain only unique coordinates for fibers from AFQ 
   def get_unique_fibers(self,CSD):

      # save CSD file name for svaing later
      fname = CSD;

      #  load in matlab AFQ style CSD bundle segmentations                                                                
      CSD = scipy.io.loadmat(CSD);

      # grab out AFQ fgs struct containing all bundle coords for 20 different bundles                                   
      fgs = CSD.get('fgs');

      # loop through each bundle coord set to grab bundle name and coords                                                
      for i in range(20):

         # grab bundle name                                                                                              
         bundle_name = str(list(fgs[0][i][0])[0]);

         # grab bundle coords array with one array per fiber                                                             
         bundle_coords = fgs[0][i][8];

         # get number of fibers for current bundle                                                                       
         num_fibers = len(bundle_coords);


         # remove repeat 3d coordinates for each fiber of each bundle                                              
         for j in range(num_fibers):

            # grab xyz coords for current fiber
            fiber_coords = bundle_coords[j][0];
            fiber_coords = np.transpose(fiber_coords); # flip to nx3

            # create array of index for each coordinate in fiber coords
            index_arr = range(len(fiber_coords));
            ncols = fiber_coords.shape[1]; # number of columns = 3
            dtype = fiber_coords.dtype.descr * ncols;
            struct = fiber_coords.view(dtype); # create struct for finding repeats using np.unique
            [unique_fiber_coords, idx] = np.unique(struct, return_index=True);
            idx = set(idx); # set of indexes of repeats in fiber_coords (rows)
            index_arr = set(index_arr); 
            repeats = np.array(list(index_arr.difference(idx))); # difference returns only the repeat idxs
            if len(repeats) !=0:
               fiber_coords[repeats] += 0.1 # add 0.1 to make sure coords are slightly different
            unique_fiber_coords = unique_fiber_coords.view(fiber_coords.dtype).reshape(-1, ncols);
            unique_fiber_coords = np.transpose(unique_fiber_coords);
            fiber_coords = np.transpose(fiber_coords);
            bundle_coords[j][0] = fiber_coords;


         #  replace old bundle coordinates with new unique only coordinates
         fgs[0][i][8] = bundle_coords;
   
      # save the new CSD file with unique fiber coordinates
      scipy.io.savemat(fname, {'fgs':fgs});


   # member function to load in results of AFQ wm bundle segmentation 
   # and shift fiber coords into fs sRAS 
   def parse_CSD_seg(self, CSD,c_RAS):


       # load in matlab AFQ style CSD bundle segmentations
       CSD = scipy.io.loadmat(CSD);

       # load in c_RAS shift vector
       c_RAS = scipy.io.loadmat(c_RAS).get('c_RAS');
    
       # grab out AFQ fgs struct containing all bundle coords for 20 different bundles
       fgs = CSD.get('fgs');
   
       # loop through each bundle coord set to grab bundle name and coords
       for i in range(20):
       
          # grab bundle name
          bundle_name = str(list(fgs[0][i][0])[0]);
        
          # grab bundle coords array with one array per fiber
          bundle_coords = fgs[0][i][8];

          # get number of fibers for current bundle
          num_fibers = len(bundle_coords);

          # apply c_RAS shift to each 3d coord for each fiber of each bundle
          for j in range(num_fibers):
       
             fiber_coords = bundle_coords[j][0];
             fiber_coords = np.transpose(fiber_coords);
             fiber_coords = fiber_coords - c_RAS;
             bundle_coords[j][0] = fiber_coords;

       
          # save current bundle in .mat file
          fname = bundle_name + '_sRAS.mat';
          scipy.io.savemat(fname, {'fg':bundle_coords});



   # member function for applying CT 2 MRI registration to
   # electrode center coordinates from 3d Hough transform
   def reg_hough_coords(self, coords_img,source, target):

      '''Member function for applying CT 2 MRI registration to
         electrode center coordinates from 3d Hough transform localization.
         
         Function takes as input:
 
         1.) coords img - a NIFTI image with 1's present
             in the position of CT electrode centers.
         2.) source - the patients native CT  NIFTI image
             used to calculate the hough gradient.
         3.) target - the patients orig.nii freesurfer volume
             output in the /mri folder after running recon-all.'''    


      # load in detected hough electrode coordinates
      # from source CT img
      coords_img = nib.load(coords_img);
      coords_data = coords_img.get_data();
    
      # grab 1's from img representing detected 
      # electrode centers from hough transform
      coords = np.nonzero(coords_data);
      x = np.array(coords[0]);
      y = np.array(coords[1]);
      z = np.array(coords[2]);
      vox_coords = np.column_stack((x,y,z));
    
      # save hough electrode centers as voxel coords in native CT space
      scipy.io.savemat(self.CT_dir + '/Hough_vox_coords.mat', {'elecmatrix':vox_coords});
    
      # load in CT2MM affine matrix
      CT_img = nib.load(source);
      CT_data = CT_img.get_data();
      CT2mm = CT_img.get_affine();
    
      # load in target MRI affine
      MRI_img = nib.load(target);
      MRI_data = MRI_img.get_data();
      MRI_hdr = MRI_img.get_header();
      MRI2mm = MRI_img.get_affine();
      CT2MRI = np.linalg.inv(MRI2mm).dot(CT2mm);
    
      # apply CT2MRI affine registration
      reg_coords = nib.affines.apply_affine(CT2MRI,vox_coords);
      scipy.io.savemat(self.CT_dir + '/Hough_reg_vox_coords.mat', {'elecmatrix':reg_coords});
    
      # save hough reg voxel coords in a NIFTI img for QC
      elecs_img = np.zeros(MRI_data.shape);
      for i in reg_coords:
         elecs_img[i[0], i[1], i[2]] = 1;
      N = nib.Nifti1Image(elecs_img,MRI2mm,MRI_hdr);
      N.to_filename(self.CT_dir + '/Hough_reg_vox_coords.nii.gz');
    
      # vox2Ras affine for all fs orig volumes
      vox2RAS = np.array([[  -1.,    0.,    0.,  128.],
            [   0.,    0.,    1., -128.],
            [   0.,   -1.,    0.,  128.],
            [   0.,    0.,    0.,    1.]]);
    
      # apply vox2RAS to get surface mesh (sRAS) coordinates      
      sRAS_coords = nib.affines.apply_affine(vox2RAS,reg_coords);
      scipy.io.savemat(self.CT_dir + '/Hough_sRAS_coords.mat', {'elecmatrix':sRAS_coords});


      return sRAS_coords;


   # function to run nmi coregistration between two NIFTI images
   # using SPMs normalizated mutual information algorithm.
   def reg_img(self,source,target):

      '''function to run nmi coregistration between two NIFTI images
         using SPMs normalizated mutual information algorithm.'''


      # create instance of spm.Corregister to run mutual information
      # registration between source and target img
      coreg = spm.Coregister();
    
      # set source and target inputs to coregister
      coreg.inputs.source = source;
      coreg.inputs.target = target;
   
      # set algorithm to normalized mutual information
      coreg.inputs.cost_function = 'nmi';

      # run coregistration
      print "::: Computing registration between CT and MRI :::";
      coreg.run();


   
   # function to run matlab project_electrodes code
   # for projecting electrode coordinates onto a surface mesh
   def project_electrodes(self, mesh, elecs, hem, out):
       
       # compute hough transform to find electrode center coordinates
       mlab = matlab.MatlabCommand();
       mlab.inputs.script = "addpath(genpath('%s')); \
                             path_to_imaging_data = '/Volumes/Thunderbolt_Duo/imaging/subjects/'; \
                             load('%s'); \
                             load('%s'); \
                             hem = %d; '\
                             debug_plots = 0; \
                             [elecs_proj] = project_electrodes(cortex, \
                             elecmatrix, hem, debug_plots); \
                             elecmatrix = elecs_proj; \
                             save('%s', elecmatrix);" %(self.project_elecs_dir, mesh, elecs, hem, out);
                             
       print "::: Loading Mesh data :::";
       print "::: Projecting electrodes to mesh :::";
       out = mlab.run();
       print "::: Done :::";
       
       return out;


   # function to use fsl FNIRT non linear registration
   # to warp a patients normalized freesurfer brain (skullstripped)
   # to freesurfer CVS MNI152 space
   def warp_FNIRT(self, source, target):


      '''function to use fsl FNIRT non linear registration                                                            
         to warp a patients normalized freesurfer brain (skullstripped)                                                        to freesurfer CVS MNI152 space'''

      
      # create FNIRT instance and set smoothing, sampling, and resolution scheme
      reg = fsl.FNIRT();
      reg.inputs.in_fwhm = [8, 4, 2, 2];
      reg.inputs.subsampling_scheme = [4, 2, 1, 1];
      reg.inputs.warp_resolution = (6, 6, 6);

      # name of outfile and coeficients file containing warping params
      out = source.replace('.nii', '_warped.nii');
      coefs = source.replace('.nii', '_coefs.nii');
    
      # begin FNIRT non linear reg
      res = reg.run(in_file=source, ref_file=target, warped_file=out, fieldcoeff_file=coefs);
