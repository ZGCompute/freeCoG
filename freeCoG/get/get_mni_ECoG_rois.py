from __future__ import division
import os
import numpy
import scipy.io
from subprocess import call

# Define mmToVox transformation
def mmToVox(mmcoords):
	#function to convert mm coordinates in the standard 2mm MNI atlas into voxel coordinates
	voxcoords = ['','','']
	voxcoords[0] = str(int(round(int(mmcoords[0])/2))*-1+45)
	voxcoords[1] = str(int(round(int(mmcoords[1])/2))+63)
	voxcoords[2] = str(int(round(int(mmcoords[2])/2))+36)
	return voxcoords
	
# setup fsl
command = 'FSLDIR=/usr/local/fsl';
print command;
call(command,shell=True);
command = '${FSLDIR}/etc/fslconf/fsl.sh';
print command;
call(command,shell=True);
command = 'PATH=${FSLDIR}/bin:${PATH}';
print command;
call(command,shell=True);
command = 'export FSLDIR PATH';
print command;
call(command,shell=True);

# Nav to dir where 
mni_dir = '/Volumes/Thunderbolt_Duo/imaging/subjects/EC77/elecs/mni/';
os.chdir(mni_dir);

# set ROI params 
radius = 2;
				
# load .mat file with all mni electrodes
fname = 'mni_depth_probes.mat';
mnicoords = scipy.io.loadmat(fname);
mnicoords = mnicoords.get('mni_elcoord');

# make a dir for putting all the rois in
os.mkdir('mni_vox_rois')
os.chdir('mni_vox_rois')

# For electrode 
i=1;
for coord in mnicoords:

    # set x,y, and z for coords 2 vox
    x = coord[0];
    y = coord[1];
    z = coord[2];

    # grab coords for 1st electrode
    coords = [x,y,z];
    voxcoords = mmToVox(coords);

    # create roi
    os.system("fslmaths $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask -roi %s 1 %s 1 %s 1 0 1 tmp" % (voxcoords[0],voxcoords[1],voxcoords[2]));
    os.system("fslmaths tmp -kernel sphere %s -fmean tmp" % (radius));
    b = str(i);
    outfile = ('Elec%s_roi.nii.gz' % (b));
    os.system("fslmaths tmp -thr .00001 -bin %s" % outfile);
    i +=1;

