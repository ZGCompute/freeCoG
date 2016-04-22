# elecs2brainVox.py
# function takes as input electrode surface-RAS coordinates, 
# an fs orig volume, and outputs voxel coordinates
# and nifti volume with the electrodes as 1's in orig voxel space

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

def elecs2brainVox(img, elecs_file, voxCoords_fname, elecsVol_fname):

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
