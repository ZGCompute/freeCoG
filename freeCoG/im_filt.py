import nibabel as nib
import nipy
from nipy.algorithms.kernel_smooth import LinearFilter
import numpy as np

def img_filt(fname):
    
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
