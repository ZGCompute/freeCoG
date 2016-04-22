# run_elecDetect.py

# script calls CT_electhresh to threshold post-implant CTs
# to increase SNR for CT voxels with platinum ECoG electrodes 
# present. It then calls im_filt.py to smooth the thresholded CT
# with a gaussian kernel before running Hough transform circle detection
# and subsequent kmeans clustering to find cluster centroids
# for electrode localization

import os

#navigate to python script repository to import CT scripts
scripts_dir = "/Volumes/Thunderbolt_Duo/Scripts/Python/";
os.chdir(scripts_dir);

from CT_electhresh import *
from im_filt import *
from run_HoughDetect import *
from Kmeans_cluster import *
from mask_CTBrain import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors

def elecDetect(CT_dir, T1_dir,img, thresh, frac, num_gridStrips, n_elects):
    
    #navigate to dir with CT img
    os.chdir(T1_dir);
    
    # Extract and apply CT brain mask
    mask_CTBrain(CT_dir,T1_dir, frac);
    masked_img = 'masked_' + img;
    
    # run CT_electhresh
    print "::: Thresholding CT at %f stdv :::" %(thresh);
    electhresh(masked_img,thresh);
    
    # run im_filt to gaussian smooth thresholded image
    fname = "thresh_" + masked_img;
    print "::: Applying guassian smooth of 2.5mm to thresh-CT :::";
    img_filt(fname);
    
    # run hough transform to detect electrodes
    new_fname = "smoothed_" + fname;
    print "::: Applying 2d Hough Transform to localize electrode clusts :::";
    elect_coords = CT_hough(CT_dir,new_fname);
    
    # set up x,y,z coords to view hough results
    x = elect_coords[:,0];
    x = x.tolist();
    y = elect_coords[:,1];
    y = y.tolist();
    z = elect_coords[:,2];
    z = z.tolist();
    
    # plot the hough results
    s = 64;
    fig = plt.figure();
    ax = fig.add_subplot(111,projection='3d');
    Axes3D.scatter3D(ax,x,y,z,s=s);
    plt.show(fig);

    # run K means to find grid and strip clusters
    n_clusts = num_gridStrips; # number of cluster to find = number of electrodes
    iters = 1000;
    init_clusts = n_clusts;
    print "::: Performing K-means cluster to find %d grid-strip and noise clusters :::" %(n_clusts);
    [clust_coords, clust_labels] = clust(elect_coords, n_clusts, iters, init_clusts);
    
    # run K means to find electrode center coordinates
    n_clusts = n_elects; # number of cluster to find = number of electrodes
    iters = 1000;
    init_clusts = n_clusts;
    print "::: Performing K-means cluster to find %d electrode centers :::" %(n_clusts);
    [center_coords, labels] = clust(elect_coords, n_clusts, iters, init_clusts);
    
    print "::: Finished :::";

    return elect_coords, clust_coords, clust_labels





