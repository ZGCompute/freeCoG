# get_elecLabels.py
# function to obtain anatomical labels of electrodes 
# based on euclidean distance from freesurfer segmented
# pial roi meshes

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: February 18, 2015

import scipy
import scipy.io
import os
from scipy.spatial.distance import cdist
import numpy as np
import collections


def get_elecLabels(subj, elecs, hem):
    
    # 1.) load in elecs and mesh data
  
    # nav to subject root dir
    root = '/Volumes/Thunderbolt_Duo/imaging/subjects/' + subj;
    os.chdir(root);
    
    # load electrode coordinates
    os.chdir(root + '/elecs');
    elec_name = elecs.strip('.mat');
    elecs = scipy.io.loadmat(elecs);
    elecs = elecs.get('elecmatrix');
    
    # load pial mesh
    rois_dir = (root + '/fs_data/surf/%s_pial_rois/vert' %(hem));
    os.chdir(rois_dir);

    # 2.) form roi list and data containers
   
    #create a list of rois to iterate over
    roi_list = [i for i in os.listdir(os.getcwd()) if (i.endswith('.mat') and i.startswith('ctx-'+ hem))];
    
    # empty dict for all labeled elecs per roi
    elecs_all = {}; 
    
    # empty dict for all euclidean distances to closest mesh point
    # for each elec (values) per roi (key)
    min_eucDists_all = {};
    
    # dict to store rois as keys, and a list of all electrode indices 
    # for electrodes close in euclidean space to that roi mesh
    roi_elec_inds_all = {};
    
    # convert coords to list. for grabbing indices of electrodes using .index
    elecs_list = list(elecs); # list with each coordinate as an XYZ array
    elecs_list = [list(i) for i in elecs_list];# list with each coordinate as an XYZ list
    
    # iterable list to iterate over each electrode in list
    elec_nums = np.linspace(0,len(elecs)-1, num=(len(elecs)));
 
    # 3.) iterate over roi list to find closest gyrul mesh of each electrode
  
    for roi in roi_list:
        
       # get vertex data for current roi
       verts = scipy.io.loadmat(roi);
       verts = verts.get('vert');
    
       # calculate euclidean distance between all electrodes and all mesh points
       eucDist = cdist(elecs, verts, 'euclidean');
    
       # get the minimum euc mesh distance for each electrode and current roi
       min_eucDists = [min(eucDist[i,]) for i in elec_nums];
       min_eucDists = scipy.array(min_eucDists);
    
       # collect all electrodes with a euc distance < 5 to current roi mesh
       roi_elecs = [elecs[i,:] for i in elec_nums if min_eucDists[i] <= (0.55 *scipy.std(min_eucDists))];
       roi_elecs = np.array(list(roi_elecs));
    
       # add selected electrodes to dict entry for current roi
       # elecs_all is a dict with rois as keys and lists of elec coords with 
       # close euclidean distances as values
       elecs_all[roi] = roi_elecs;
       
       # get indices of all electrodes selected for current roi       
       # empty list for poulating with electrode inds
       roi_elec_list = list(roi_elecs);
       roi_elec_list = [list(i) for i in roi_elec_list]; # convert list of arrays to list of lists[xys]
       roi_elec_inds = []; 
       elec_ind_nums = np.linspace(0,len(roi_elecs)-1, num=(len(roi_elecs)));
       elec_ind_nums = elec_ind_nums.astype('int'); # index into roi_elec_list
       
       # populate elec_ind_nums with indeces of all ellectrodes for current roi
       for i in elec_ind_nums:
           roi_elec_inds.append(elecs_list.index(roi_elec_list[i])); # index into original elecs_lsit
       
       # add roi elec inds to dict of all roi elec inds: each roi(key) has a list
       # of indices(val)
       roi_elec_inds_all[roi] = roi_elec_inds;
       
       # append euclidian distances to current mesh for all electrodes
       # min_euc_Dists_all is a dict with rois for keys, each with a list of
       # euclidean distances from each electrode (for the closest point on the
       # current mesh) 
       min_eucDists_all[roi] = min_eucDists;
       

    # 4.) find electrodes that are close to two rois and compare to find closest 
    # roi for a given electrode 
  
    # closest euclidean distance for an electrode
    # get dict of electrode indices (w), with values being a list of rois overlapping
    w = collections.defaultdict(list); 
    for k,v, in roi_elec_inds_all.iteritems():
       for i in v: w[i].append(k);
    
    # iterate over electrode indexes
    for k in w.keys():
           
       # if for a given electrode index, there are two or more rois listed with the same electrode
       if len(w[k]) >1:
          min_dist = 1000000; # set initial min dist to 1 mil
             
          # find which region has a smaller euc distance to that electrode
          for region in w[k]:
                 
             # if the current region has smaller eucdist make it the min
             if min_eucDists_all[region][k] < min_dist:
                min_dist = min_eucDists_all[region][k];
                closest = region;
                     
             # change the key label to the closest region
             w[k] = closest;
                 

    # output final list of each electrodes label
    elec_labels =[];
    for k in w.keys():
        if type(w[k]) == list:
           elec_labels.append(w[k][0]);
        else:
           elec_labels.append(w[k]);
            
    # get rid of 'ctx-' and '_vert.mat' for each label
    elec_labels = [label.strip('_vert.mat') for label in elec_labels];
    elec_labels = [label.strip('ctx-') for label in elec_labels];
    
    # list of all rois included in final labeling
    label_rois = set(elec_labels);
    label_rois = list(label_rois);
        
    
    # 5.) Save final outputs
  
    # create dirs for saved labeled elecs
    os.chdir(root + '/elecs');
    elecs_dir = os.getcwd();
    labels_dir = elecs_dir + '/labels';
    
    # check if labels dir already exists
    if os.path.isdir(labels_dir) == True:
        
        None;
        
    else:
        
       os.mkdir('labels');
    
    # go in to labels dir   
    os.chdir(labels_dir);
    
    # save coords .mat for each roi
    elec_nums = [int(num) for num in elec_nums];
    
    #check if any electrodes went missing and didn't get labeled
    elec_nums_set = set(elec_nums);
    missing = elec_nums_set.difference(w.keys());
    missing = list(missing);
    
    # if there are missing electrodes add the missing electrode indeces to
    # the list of electrode keys, and fill in the label using the adjacent 
    # electrode labels in elec_labels
    if missing != None:
       keys = w.keys()
       for missed in missing:
           keys.insert(missed, missed);
           elec_labels.insert(missed,elec_labels[missed]);
           
    
       keys.sort(); # make sure keys is in order
          
           
        
    # save each list of labeld elecs
    for roi in label_rois:
        
        # container for elecs on current roi to save
        labeled_elecs = [];
        
        # add all electrodes into container if label is current roi
        for i in elec_nums:
            if elec_labels[i] == roi:
                labeled_elecs.append(elecs[i]);
        
        # save a .mat for the current roi coords
        labeled_elecs = np.array(labeled_elecs);
        fname = elec_name + '_' + roi + '_elecs.mat';
        scipy.io.savemat(fname, {'elecmatrix':labeled_elecs});


    # return list of elec_labels
    return elec_labels;
    
            
                   
 
       
       
       
       
    
    
    
    
    