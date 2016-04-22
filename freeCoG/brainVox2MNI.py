# brainVox2MNI.py
# function takes as input electrode voxel coordinates in 
# freesurfer orig space, and a non-linear transformation
# matrix for MNI warping, and outputs voxel and freeusrfer
# Surface-RAS coordinates in .nii volume and .mat 3d coords

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: January 30, 2015


import nibabel as nib
import scipy.io
import numpy
import os
import subprocess
import itertools

def brainVox2MNI(elec_vox, norm_mat, MNI_vol, MNI_RASmat_fname, MNI_voxMat_fname, MNI_elecsVolFname):
    
   # load transform into numpy array
   norm_mat = scipy.io.loadmat(norm_mat);
   norm_mat = norm_mat.get('Affine');
   
   # load elec voxel coords
   elec_vox = scipy.io.loadmat(elec_vox);
   elec_vox = elec_vox.get('elecmatrix');
   
   # prep intercept for applying the warp
   intercept = numpy.ones(len(elec_vox));
   mni_elecs  = numpy.column_stack((elec_vox,intercept));
   mni_elecs  = numpy.transpose(mni_elecs);
   mni_elecs  = numpy.matrix.dot(norm_mat,mni_elecs);

   
   # apply vox2RAS affine to get fs surface coords
   affine = numpy.array([[  -1.,    0.,    0.,  128.],
          [   0.,    0.,    1., -128.],
          [   0.,   -1.,    0.,  128.],
          [   0.,    0.,    0.,    1.]]);
   mni_RAS = numpy.matrix.dot(affine, mni_elecs);
   mni_RAS = numpy.transpose(mni_RAS);
   mni_RAS = mni_RAS[:,0:3];
   
   # save mat for mni surface RAS coords
   scipy.io.savemat(MNI_RASmat_fname, {'elecmatrix':mni_RAS});
   
   # save mat with mni voxel coordinates
   mni_elecs = numpy.transpose(mni_elecs);
   mni_elecs =mni_elecs[:,0:3];
   mni_elecs = mni_elecs.astype('int');
   scipy.io.savemat(MNI_voxMat_fname, {'elecmatrix':mni_elecs});
   
   # save nifti volume in MNI space with electrodes as 1s
   MNI_brain = nib.load(MNI_vol);
   dat = MNI_brain.get_data();
   
   # empty img for placing elecs in
   MNI_elecs_img = numpy.zeros(dat.shape);
   
   # put ones in the vox coords of elecs in empty img
   for i in mni_elecs:
      MNI_elecs_img[i[0], i[1], i[2]] = 1;
      
   # write mni vox image to nifti
   hdr = MNI_brain.get_header();
   affine = MNI_brain.get_affine();
   N = nib.Nifti1Image(MNI_elecs_img,affine,hdr);
   N.to_filename(MNI_elecsVolFname);
   
   return mni_RAS, mni_elecs;
   
 
