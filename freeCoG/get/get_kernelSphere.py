# get_CT_kernel_sphere.py
# function takes as input a thresholded CT image
# and ouputs a spherical kernel image of the same size

import nibabel
import os

def get_kernel(img):
    
   # init vox coords for kernel center
   voxcoords = [100, 100, 100];

   # load in image data
   im = nibabel.load(img);
   im_dat = im.get_data();
   
   # define kernel radius
   radius = 3; # 2 voxel radius
   
   # use fslmaths to create sphere of 1's at voxcoord
   outfile  = 'kernel_sphere.nii.gz';
   print "::: Creating Kernel Sphere with fslmaths :::";
   os.system("fslmaths %s -mul 0 -add 1 -roi %s 1 %s 1 %s 1 0 1 tmp" % (img,voxcoords[0],voxcoords[1],voxcoords[2]));
   os.system("fslmaths tmp -kernel sphere %s -fmean tmp" % (radius));
   print "::: Binarizing volumetric Kernel Sphere :::";
   os.system("fslmaths tmp -thr .00001 -bin %s" % outfile);
   
   # load in kernel to return mat
   kernel = nibabel.load(outfile);
   kernel= kernel.get_data();
   kernel = kernel[90:110,90:110,90:110];
   
   return kernel;
   
   
   
   