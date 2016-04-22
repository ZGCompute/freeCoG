# subcortFs2mlab(subcort)

# function reads in freesurfer subcortical surfaces to convert them 
# to ASCII to create matlab meshes


import os
import nipype
import scipy
import scipy.io
import numpy as np

def subcortFs2mlab(subcort,nuc):
    
    # setup freesurfer
    #os.system('export FREESURFER_HOME=/usr/local/freesurfer-5.3.0');
    #os.system('source $FREESURFER_HOME/SetUpFreeSurfer.csh'
    
    # use freesurfer mris_convert to get ascii subcortical surface
    subcort_ascii = subcort;
    #os.system('mris_convert %s %s' %(rh, rh_ascii));
    
    # clean up ascii file and extract matrix dimensions from header
    subcort = open(subcort_ascii, 'r');
    subcort_mat = subcort.readlines();
    subcort.close();
    subcort_mat.pop(0); # get rid of comments in header
    subcort_mat = [item.strip('\n') for item in subcort_mat];# get rid of new line char
    
    # extraxt inds for vert and tri    
    subcort_inds = subcort_mat.pop(0); 
    subcort_inds = subcort_inds.split(' '); # seperate inds into two strings
    subcort_inds = [int(i) for i in subcort_inds];# convert string to ints
    

    # get rows for vertices only, strip 0 column, and split into seperate strings
    subcort_vert = [item.strip(' 0') for item in subcort_mat[:subcort_inds[0]]];
    subcort_vert = [item.split('  ') for item in subcort_vert]; # seperate strings
    
    # create containers for each column in vert matrix
    x = [];
    y = [];
    z = [];
    
    # fill the containers with float values of the strings in each column
    for i in subcort_vert:
        x.append(float(i[0]));
        y.append(float(i[1]));
        z.append(float(i[2]));
        
    # convert to scipy mat
    x = scipy.mat(x);
    y = scipy.mat(y);
    z = scipy.mat(z);
    
    # concat columns to one n x 3 matrix
    x = x.transpose();
    y = y.transpose();
    z = z.transpose();
    subcort_vert = scipy.column_stack((x,y,z));
    scipy.io.savemat('subcort_%s_vert.mat' %(nuc), {'vert':subcort_vert});# save vert matrix
    
    # get rows for triangles only, strip 0 column, and split into seperate strings
    subcort_tri = [item[:-2] for item in subcort_mat[subcort_inds[0] +1:]];
    subcort_tri = [item.split(' ') for item in subcort_tri]; # seperate strings  
    
    # create containers for each column in vert matrix
    x = [];
    y = [];
    z = [];
    
    # fill the containers with float values of the strings in each column
    for i in subcort_tri:
        x.append(int(i[0]));
        y.append(int(i[1]));
        z.append(int(i[2]));
        
    # convert to scipy mat
    x = scipy.mat(x);
    y = scipy.mat(y);
    z = scipy.mat(z);
    
    # concat columns to one n x 3 matrix
    x = x.transpose();
    y = y.transpose();
    z = z.transpose();
    subcort_tri = scipy.column_stack((x,y,z));
    subcort_tri = subcort_tri + 1;
    scipy.io.savemat('subcort_%s_tri.mat' %(nuc), {'tri':subcort_tri});# save tri matrix  

    # convert inds to scipy mat
    subcort_inds = scipy.mat(subcort_inds);
    scipy.io.savemat('subcort__%s_inds.mat' %(nuc), {'inds':subcort_inds});# save inds
    
    
    

