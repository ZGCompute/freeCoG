# get_cvsWarp.py
# Script to to compute non-linear warping of patient mri and surface rendereing
# into freesurfer cvs152_inMNI template space using mri_cvs_register
# Then uses applyMorph to apply the warping to  list of voxel coordinates
# representing electrode positions localized from post-op CT


# Author: Zachary Greenberg
# Date: 2/25/2015
# usage: import get_cvsWarp; get_cvsWarp(subj,hem)

import os
import numpy as np
import scipy 
import scipy.io
from elecs2brainVox import *
import nibabel as nib
import subprocess
import shlex

def get_cvsWarp(subj,hem):


   # define patient root and fs_dir
   fs_subjects_dir = '/crtx1/zgreenberg/subjects/';
   patient_root = fs_subjects_dir + subj;
   imaging_dir = '/data_store2/imaging/subjects/' + subj + '/';   

   #nav to patient dir with mri for warping to cvs
   os.chdir(patient_root);
   
   #1.)
   #run cvs register
   orig = patient_root + '/mri/brain.nii'; # orig is mri in fs orig space
   os.system('mri_convert %s/mri/brain.mgz %s' %(patient_root,orig)); # create brain.nii
   template = 'cvs_avg35_inMNI152'; # template to register to
   command = 'mri_cvs_register --mov %s --template %s --m3d --openmp 8 --keepelreg --hemi %s' %(subj,template,hem);
   args = shlex.split(command);
   p = subprocess.Popen(args);
   p.wait(); # wait to continue until mri_cvs finsihes 
   
   #2.)
   # apply non-linear warp computed from mri_cvs_register
   # to all vox coord lists to get warped vox and surface coords in MNI
   elecs_dir = imaging_dir + 'elecs/';
   cvs_dir = patient_root + '/cvs/';
   template_brain = fs_subjects_dir + template + '/mri/brain.mgz';
   transform = cvs_dir + 'combined_tocvs_avg35_inMNI152_elreg_afteraseg-norm.tm3d';
   mni_dir = elecs_dir + 'cvs_mni/';
   os.mkdir(mni_dir);
   
   #loop through all patients electrode files to warp each into MNI152
   for coord_list in os.listdir(elecs_dir):
   
      if coord_list.endswith('.mat') and coord_list.find('interp') == -1 and coord_list.find('corners') == -1 and coord_list.find('vox') == -1:

         # load in .mat file with RAS coords and convert to voxel coordinates and nifti volume
         voxCoords_fname = elecs_dir + (coord_list.strip('.mat')) + '_vox.mat';
         elecsVol_fname = voxCoords_fname.replace('.mat', '.nii.gz');
         elecs_file = elecs_dir + coord_list;
         elecs2brainVox(orig, elecs_file, voxCoords_fname, elecsVol_fname);

         # load voxel coords .mat and save as .txt. for applyMorph
         voxel_coords = scipy.io.loadmat(voxCoords_fname);
         voxel_coords = voxel_coords.get('elecmatrix');
         voxCoords_txt = voxCoords_fname.replace('.mat', '.txt');
         np.savetxt(voxCoords_txt, voxel_coords, delimiter=' ', fmt='%-7.1f');

         # apply transform
         mni_voxCoords = mni_dir + 'mni_' + coord_list.replace('.mat', '_vox.txt'); # outfile
         os.system('applyMorph --template %s --transform %s point_list %s %s linear' %(template_brain,transform,voxCoords_txt, mni_voxCoords));
         
      
         # save each as .mat for visualization
         coords = np.loadtxt(mni_voxCoords);
         coords = np.array(coords);
         coords = coords[0:len(coords)-1];
         brain_incvs = cvs_dir + 'final_CVSmorphed_tocvs_avg35_inMNI152_norm.mgz';
         cvs_brain = nib.load(brain_incvs);
         dat = cvs_brain.get_data();
         
         # empty img for placing elecs in
         MNI_elecs_img = np.zeros(dat.shape);
         
         # put ones in the vox coords of elecs in empty img
         for i in coords:
            MNI_elecs_img[i[0], i[1], i[2]] = 1;

         # write mni vox image to nifti
         hdr = cvs_brain.get_header();
         affine = cvs_brain.get_affine();
         N = nib.Nifti1Image(MNI_elecs_img,affine,hdr);
         MNI_elecsVolFname = mni_voxCoords.replace('.txt', '.nii.gz');
         N.to_filename(MNI_elecsVolFname);

         # save .mat with vox coords
         new_fname = mni_voxCoords.replace('.txt', '.mat');
         scipy.io.savemat(new_fname, {'elecmatrix': coords});

         # save as nifti volume and RAS coords for visualizing mni result
         # prep intercept for applying the ras transform
         intercept = np.ones(len(coords));
         mni_elecs  = np.column_stack((coords,intercept));
         mni_elecs  = np.transpose(mni_elecs);
         
         # apply vox2RAS affine to get fs surface coords
         affine = np.array([[  -1.,    0.,    0.,  128.],
             [   0.,    0.,    1., -128.],
             [   0.,   -1.,    0.,  128.],
             [   0.,    0.,    0.,    1.]]);  
     
         mni_RAS = np.matrix.dot(affine, mni_elecs);
         mni_RAS = np.transpose(mni_RAS);
         mni_RAS = mni_RAS[:,0:3];
         
         # save mat for mni surface RAS coords
         MNI_RASmat_fname = mni_voxCoords.replace('.txt', '.mat');
         scipy.io.savemat(MNI_RASmat_fname, {'elecmatrix':mni_RAS});
   
   # clean up and organize outputs
   os.chdir(elecs_dir);
   os.system('mkdir nii/');
   os.system('mv *.nii* nii/');
   os.system('mkdir vox/');
   os.system('mv *vox.* vox/');
   os.chdir(mni_dir);
   os.system('mkdir nii/');
   os.system('mv *.nii* nii/');
   os.system('mkdir vox/');
   os.system('mv *vox.* vox/');
       



