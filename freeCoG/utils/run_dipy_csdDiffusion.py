# run_dipy_csdDiffusion.py

# function to run constrained spherical deconvolution fiber tracking 
# on HARDI diffusion weighted imaging

import os
import dipy
import nibabel as nib
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel as CsdModel
from dipy.reconst.csdeconv import auto_response
from dipy.data import get_sphere
import numpy as np
#from reconst_csa import *
from dipy.reconst.interpolate import NearestNeighborInterpolator
from dipy.tracking.markov import (BoundaryStepper,FixedSizeStepper,ProbabilisticOdfWeightedTracker)
from dipy.tracking.utils import seeds_from_mask
from dipy.segment.mask import median_otsu
import matplotlib.pyplot as plt
from dipy.viz import fvtk
from dipy.viz.colormap import line_colors

# init subjects id
subj = 'EC71';
os.chdir('/Users/Zach/Desktop/Data/MR_data/%s/DTI_55_2000/raw/' %(subj));

# load diffusion data, bvals, and vecs
img = nib.load('dwi.nii')
data = img.get_data();
affine = img.get_affine();
header = img.get_header();
voxel_size = header.get_zooms()[:3];

# slice out b0
b0  = data[:,:,:,0];
img2 = nib.Nifti1Image(b0,affine);
nib.save(img2, "b0.nii.gz");

# load bvals bvecs
fbval = "dwi.bval";
fbvec = "dwi.bvec";
bvals, bvecs = read_bvals_bvecs(fbval, fbvec);
gtab = gradient_table(bvals, bvecs);

# fit constrained spherical deconvolution model
response, ratio = auto_response(gtab, data, roi_radius=10, fa_thr=0.7);
csd_model = CsdModel(gtab, response);
csd_fit = csd_model.fit(data);

# init fiber odf
sphere = get_sphere('symmetric724');
csd_odf = csd_fit.odf(sphere);
stepper = FixedSizeStepper(1);


# load in b0 and create mask
b0_img = nib.load('b0.nii.gz');
b0_data = b0_img.get_data();
b0_mask, mask = median_otsu(b0_data, 2, 1);

# save mask for b0
mask_img = nib.Nifti1Image(mask.astype(np.float32), b0_img.get_affine())
b0_img = nib.Nifti1Image(b0_mask.astype(np.float32), b0_img.get_affine())
fname = 'b0';
nib.save(mask_img, fname + '_binary_mask.nii.gz');
nib.save(b0_img, fname + '_mask.nii.gz');

# init seeds from mask of b0
zooms = header.get_zooms()[:3];
seeds = seeds_from_mask(mask, [1, 1, 1], zooms);
seeds = seeds[:2000];
maskdata = mask.data;

# run streamline probobalistic fiber tracking
interpolator = NearestNeighborInterpolator(maskdata, zooms);
pwt = ProbabilisticOdfWeightedTracker(csd_model, interpolator, mask,stepper, 20, seeds, sphere);
csa_streamlines = list(pwt);

# save streamlines
hdr = nib.trackvis.empty_header();
hdr['voxel_size'] = (2., 2., 2.);
hdr['voxel_order'] = 'LAS';
hdr['dim'] = csapeaks.gfa.shape[:3];
csa_streamlines_trk = ((sl, None, None) for sl in csa_streamlines);
csa_sl_fname = 'csa_prob_streamline.trk';
nib.trackvis.write(csa_sl_fname, csa_streamlines_trk, hdr);


# render fibers
r = fvtk.ren()
fvtk.add(r, fvtk.line(csa_streamlines, line_colors(csa_streamlines)))
