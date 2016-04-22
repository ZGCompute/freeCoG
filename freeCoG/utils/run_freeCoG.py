# run_freeCoG.py

import freeCoG
import os

# creat freeCoG instance for current patient
fcg = freeCoG();

# set MATLAB toolbox paths
fcg.VISTA_dir = '/crtx1/zgreenberg/scripts/MATLAB/vistasoft-master/';
fcg.SPM_dir = '/crtx1/zgreenberg/scripts/MATLAB/spm12/';
fcg.AFQ_dir = '/crtx1/zgreenberg/scripts/MATLAB/AFQ/';
fcg.fsLUT_dir = '/usr/local/freesurfer/';
fcg.bash_dir = '/crtx1/zgreenberg/scripts/bash/';

# run structual processing
fcg.subj = 'EC86';
fcg.hem = 'lh';
fcg.prep_recon();
fcg.get_recon();
fcg.convert_fs2mlab();
fcg.get_subcort();
fcg.get_cortSegMeshes();
fcg.get_cRAS();
fcg.make_CortROIs();

# run CT electrode detection

img = fcg.CT_dir + '/CT.nii';
fcg.electhresh(img,thresh_f=3,thresh_c=4)
fcg.img_filt()
fcg.Spherical_Hough();
coords_img = fcg.CT_dir + '/CT_elecCenters.nii';
source = fcg.CT_dir + '/CT.nii';
target = fcg.fs_dir + '/mri/orig.nii';
fcg.reg_img(source,target);
fcg.reg_hough_coords(coords_img, source, target);
fcg.clean_Hough3D(); # run locally to segment electrode arrays manually
fcg.get_clinGrid();

# loop through all patients electrode arrays to compute anatomical labels
elecs_list = os.listdir(fcg.elecs_dir);
for elecs in elecs_list:
    fcg.get_elecLabels(elecs);
    
# get MNI coords
fcg.get_cvsWarp();

# run DWI tracography
fcg.DWI_bval = 2000;
fcg.prep_DWI_brain();
fcg.DWI_preproc();
fcg.get_DWImask();
fcg.get_mrtrix_metrics();
CSD = '/path/to/CSD/seg.mat';
c_RAS = '/path/to/c_RAS.mat';
fcg.get_unique_fibers(CSD);
fcg.clean_fiberOutliers(CSD);
fcg.parse_CSD_seg(CSD,c_RAS)
fcg.get_Connectome();
