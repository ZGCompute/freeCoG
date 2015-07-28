# Get_Vertex_thickness.sh
# Description: This script will run vertex-wise thickness comparisons between two groups 
# using Freesurfer v 5.X. Subjects should all have been processed with recon-all. You will
# want to replace <variable_name> with the variable you are using to seperate your two groups
# in the code below (i.e., pathology). 

# The next step is to create .fsgd (Freesurfer Group Discriptor) file that contains information 
# about the discrete and continuous variables you want to model using linear regression.  
# This is a simple text file, where the first line is “GroupDescriptorFile1”. 
# The next two rows have two Class entries (your two groups). For example 
# if you want to compare patients with Primary Progressive Aphasia (PPA), you could enter CLASS PPA 
# for the second row, and CLASS Control for the third row of the file. The 4th row will be “Variables” 
# in this example, and will be fallowed by a row for each subject, prefixed with “Input”, fallowed by name 
# of subject in the directory, class (or group, i.e., ‘Control’), and any continuous variables you want to 
# enter into the model. See http://surfer.nmr.mgh.harvard.edu/fswiki/Fsgdf2G0V for help.

# Author: Zachary Greenberg
# UCSF Imaging CORE
# Email: Zgreenberg@memory.ucsf.edu
# Last Edited: 3/9/2014

#!/bin/bash

# Set up Environment for freesurfer use
iam=`whoami` #user name

# establish where your subject data is. SUBJECTS_DIR and FREESURFER_HOME
# will need to be changed for different workstations
export SUBJECTS_DIR=~/Applications/freeesurfer/subjects/<folder where your subjects are >
export FREESURFER_HOME=~/Applications/freeesurfer/
source FreeSurferEnv.sh

# navigate to freesurfer subjects dir
cd $SUBJECTS_DIR

# write the inital glm perameters. substitute a group variable-name for group. 
# You will want to change the contrasts based on your hypothesis and variables entered
# into the model. See http://surfer.nmr.mgh.harvard.edu/fswiki/FsgdExamples for help
python "f = open('lh-Avg-Thickness-group-cor.mtx', 'wt', encoding='utf-8'); s = '+0.0 +0.0 +0.5 +0.5'; f.write(s);")

# run mris_preproc to:
# 1. Resample each subject's data into a common space, and
# 2. Concatenate all the subjects' into a single file
# 3. Spatial smoothing (can be done between 1 and 2)

export FWHM_FILE=$SUBJECTS_DIR/glm/lh.<variable_name>.thickness.10.mgh
if [ -d "$FWHM_FILE" ]; then 
	mris_preproc --fsgd <variable_name>.fsgd \
	--cache-in thickness.fwhm10.fsaverage \
	--target fsaverage --hemi lh \
	--out lh.<variable_name>.thickness.10.mgh
	
# if data haven't been cached
else
	mris_preproc --fsgd <variable_name>.fsgd \
	--cache-in thickness.fwhm10.fsaverage \
	--target fsaverage --hemi lh \
	--out lh.<variable_name>.thickness.10.mgh	
fi

# run GLM analysis 
mri_glmfit \ 
--y lh.<variable_name>.thickness.10.mgh \ 
--fsgd <variable_name>.fsgd dods\ 
--C lh-Avg-thickness-<variable_name>-Cor.mtx \ 
--surf fsaverage lh \ 
--cortex \ 
--glmdir lh.<variable_name>.glmdir

# View the uncorrected significance map with tksurfer:
tksurfer fsaverage lh inflated \ 
-annot aparc.annot -fthresh 2 \ 
-overlay lh.<variable_name>.glmdir/lh-Avg-thickness-<variable_name>-Cor/sig.mgh

# Run monte-carlo sim to account for multiple comparisons correction
mri_glmfit-sim \ 
--glmdir lh.<variable_name>.glmdir \ 
--sim mc-z 5 4 mc-z.negative \ 
--sim-sign neg \ 
--overwrite

# plot the corrected overlay
tksurfer fsaverage lh inflated \ 
-annot lh.<variable_name>.glmdir/lh-Avg-thickness-<variable_name>-Cor/mc-z.neg4.sig.ocn.annot \ 
-fthresh 2 -curv -gray
