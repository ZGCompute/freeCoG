# loadSave_hough_dat.py 
# script to save spherical hough data as volume and get 3d coords

import scipy.io
import nibabel as nib
import nipy
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors


# load pre-processed CT to get hdr and affine
fname = 'smoothed_thresh_CT.nii'
img = nib.load(fname);
hdr = img.get_header();
affine = img.get_affine();

# load 3d hough centers and spheres
img_data = scipy.io.loadmat('elec_spheres.mat');
img_data = img_data.get('sphere_img');
new_fname = 'CT_elecSpheres.nii';
N = nib.Nifti1Image(img_data,affine,hdr);
N.to_filename(new_fname);

# repeat for center img data
img_data = scipy.io.loadmat('elec_centers.mat')
img_data = img_data.get('center_img')
new_fname = 'CT_elecCenters.nii'
N = nib.Nifti1Image(img_data,affine,hdr);
N.to_filename(new_fname);

# get img dims
z_len = img_data[0,0,:].shape;
x_len = img_data[:,0,0].shape;
y_len = img_data[0,:,0].shape;

# get coords of voxels containing center data
coords = np.nonzero(img_data);
x = coords[0];
y = coords[1];
z = coords[2];

# concatonate xyz vectors into 3 x n matrix
stackem = (x,y,z);
coords = np.vstack(stackem);
coords = np.transpose(coords);
scipy.io.savemat('3dHough_coords.mat', {'coords':coords});


# plot the hough center results
s = 202
fig = plt.figure(facecolor="white");
ax = fig.add_subplot(111,projection='3d');
Axes3D.scatter3D(ax,x,y,z,s=s);
plt.show(fig);





