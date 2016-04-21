# get_pial_roiMesh(subj, hem)

# function converts .dpv aparc from Freesurfer to matlab meshes of each
# roi on the pial surface


from pialROI_fs2mlab import *
import os

def get_pial_roiMesh(subj, hem):
   
   # location of freesurfer look up table for naming ROIs
   fsLUT_dir= '/usr/local/freesurfer/';
   os.chdir(fsLUT_dir);
    
   # open the LUT text file and grab list of label #'s and names
   fname = 'FreeSurferColorLUT.txt';
   file = open(fname, 'r');
   LUT = file.readlines();
   LUT = LUT[429:502]; # grab only crtx ROIs from .txt LUT
   LUT = [row.split("    ") for row in LUT];
   LUT = [row[0:2] for row in LUT]; # only retain label # and name
   LUT.sort();
   LUT.pop(0); # get rid of non-sense string 
   LUT.pop(0); # and get rid of unknown roi
   LUT.pop(35); # rh unknown
      
   # nav to dir where *hem.pial is located
   os.system('cd /Users/Zach/Desktop/Data/MR_data/%s/fs_data/label' %(subj));
   
   # copy areal toolbox bash scripts to cwd
   os.system('cp /Users/Zach/Desktop/Scripts/Bash/areal/bin/* .');
   os.system('cp /Users/Zach/Desktop/Scripts/Bash/areal/share/* .'); 
   
   # convert annot.apparc to .dpv
   os.system('./annot2dpv %s.aparc.annot %s.aparc.annot.dpv' %(hem,hem));
   
   # nav to dir where *hem.pial is located
   os.system('cd /Users/Zach/Desktop/Data/MR_data/%s/fs_data/surf' %(subj));
   
   # copy areal toolbox bash scripts to cwd
   os.system('cp /Users/Zach/Desktop/Scripts/Bash/areal/bin/* .');
   os.system('cp /Users/Zach/Desktop/Scripts/Bash/areal/share/* .');
   
   # split aparc into roi .srf files for each roi
   os.system('./splitsrf %s.pial.asc /Users/Zach/Desktop/Data/MR_data/%s/fs_data/label/%s.aparc.annot.dpv %s.pial_roi' %(hem, subj, hem, hem));

   # get list of all the .srf roi files created from splitsrf
   srf_list = [i for i in os.listdir(os.getcwd()) if i.endswith('.srf')];

   # create corresponding list of 
   ascii_list = [i.replace('.srf', '.asc') for i in srf_list];
   
   # convert ascii to matlab style tri and vert for meshes
   for fname in ascii_list:
      fname = fname;
      roi = fname[12:16];
      pialROI_2mlab(fname,roi);
      
   # mv all matlab files to roi folder
   os.mkdir('pial_rois');
   os.system('mv *.mat pial_rois');
   
   # move over ascii files
   os.system('mv %s.pial_roi* pial_rois' %(hem));
   os.system('cd pial_rois');
   os.mkdir('ascii');
   os.system('mv %s.pial_roi*.asc ascii' %(hem));
   os.mkdir('inds');
   os.system('mv *inds.mat inds');
