# CT_electhresh.py
# Function thresholds a CT to pull out electrodes and remove skull/fillings

# Author: Zachary I. Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: March 19, 2015


import numpy as np
from numpy import zeros
import nibabel as nib
import scipy
import mahotas


def electhresh(img,thresh_f,thresh_c):

    
   # function thresholds a CT to pinpoint electrode locations
   img_fname = img;
   img = nib.load(img);
   img_dat = img.get_data();
   img_arr = np.array(img_dat);
   new_fname = 'thresh_' + img_fname;
   img.set_filename = new_fname;
    
   # Create a container for the thresholded CT values
   CT_thresh = zeros(img_arr.shape);
    
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

                
                
