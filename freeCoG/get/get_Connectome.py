# get_Connectome.py

# function wraps mrtrix tckgen for pair wise roi based tractography
# using freesurfer rois and constrained spherical convolution

import os
import subprocess
import shlex

def get_Connectome(subj):
    

   # set up var for DWI dir and roi dir
   DWI_dir = '/Volumes/Thunderbolt_Duo/imaging/subjects/' + subj + '/DTI_HARDI_2000_55/';
   roi_dir = '/Volumes/Thunderbolt_Duo/imaging/subjects/' + subj + '/fs_data/mri/crtx_subcort_ROIs/';
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
                      %smask.mif -nthreads 12 -maxnum 10000' %(DWI_dir, roi, roi_s, roi1, roi2, DWI_dir);
           args = shlex.split(command);
           subprocess.Popen(args);
    
   exit# clean up results              
   os.system('mv *.tck roi_roi_tcks/');
