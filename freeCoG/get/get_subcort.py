# get_subcort.py

# function tesselates and converts all subjects freesurfer subcortical 
# segmentations to matlab tri-vert format for plottings as Matlab mesh

import os
from subcortFs2mlab import *

def get_subcort(subj):
    
    bash_scripts_dir = '/crtx1/zgreenberg/scripts/bash';
    python_scripts_dir = '/crtx1/zgreenberg/scripts/python';
    subjAscii_dir = ('/crtx1/zgreenberg/subjects/%s/ascii/' %(subj));
    
    # nav to scripts dir to call aseg2srf
    os.chdir(bash_scripts_dir);
    
    # tessellate all subjects freesurfer subcortical segmentations
    os.system('./aseg2srf.sh -s "%s" -l "4 5 10 11 12 13 17 18 26 28 43 44  49 50 51 52 53 54 58 60 14 15 16"' %(subj));
    
    # nav to dir with ascii outputs for each subcortical segmented surface
    os.chdir(subjAscii_dir);
    
    # get list of all .srf files and change fname to .asc
    srf_list = [fname for fname in os.listdir(os.getcwd())];
    asc_list = [fname.replace('.srf', '.asc') for fname in srf_list];
    asc_list.sort();
    for fname in srf_list:
       new_fname = fname.replace('.srf', '.asc');
       os.system("mv %s %s" %(fname, new_fname));
    
    # convert all ascii subcortical meshes to matlab vert, tri coords
    print "::: Converting all ascii segmentations to matlab tri-vert :::";
    subcort = 'aseg_058.asc';
    nuc = 'rAcumb';
    subcortFs2mlab(subcort,nuc);
 
    subcort = 'aseg_054.asc';
    nuc = 'rAmgd';
    subcortFs2mlab(subcort,nuc);

    subcort = 'aseg_050.asc';
    nuc = 'rCaud';
    subcortFs2mlab(subcort,nuc);

    nuc = 'rGP';
    subcort = 'aseg_052.asc';
    subcortFs2mlab(subcort,nuc);

    nuc = 'rHipp';
    subcort = 'aseg_053.asc';
    subcortFs2mlab(subcort,nuc);

    subcort = 'aseg_051.asc';
    nuc = 'rPut';
    subcortFs2mlab(subcort,nuc);

    subcort = 'aseg_049.asc';
    nuc = 'rThal';
    subcortFs2mlab(subcort,nuc);

    nuc = 'rLatVent';
    subcort = 'aseg_043.asc';
    subcortFs2mlab(subcort,nuc);

    subcort = 'aseg_044.asc';
    nuc = 'rInfLatVent';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_060.asc';
    nuc = 'rVentDienceph';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_004.asc';
    nuc = 'lLatVent';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_005.asc';
    nuc = 'lInfLatVent';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_010.asc';
    nuc = 'lThal';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_011.asc';
    nuc = 'lCaud';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_012.asc';
    nuc = 'lPut';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_013.asc';
    nuc = 'lGP';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_017.asc';
    nuc = 'lHipp';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_018.asc';
    nuc = 'lAmgd';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_026.asc';
    nuc = 'lAcumb';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_028.asc';
    nuc = 'lVentDienceph';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_014.asc';
    nuc = 'lThirdVent';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_015.asc';
    nuc = 'lFourthVent';
    subcortFs2mlab(subcort,nuc);
    
    subcort = 'aseg_016.asc';
    nuc = 'lBrainStem';
    subcortFs2mlab(subcort,nuc);
       
    
    
    
    
    
    
    
    
    
    
    
    