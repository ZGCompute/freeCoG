import numpy as np
import scipy
import scipy.io
import sklearn
from sklearn import cluster
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors

def clust(elect_coords, n_clusts, iters, init_clusts):
    
   # Load resultant coordinates from Hough circles transform
   #coords = scipy.io.loadmat(elect_coords);
   #dat = coords.get('elect_coords');
   dat = elect_coords;

   # Configure Kmeans
   cluster = sklearn.cluster.KMeans();
   cluster.n_clusters= n_clusts;
   cluster.init = 'k-means++';
   cluster.max_iter = iters;
   cluster.verbose = 0;
   cluster.n_init = init_clusts;
   cluster.fit(dat);

   # Grab vector for plotting each dimension
   x = list(cluster.cluster_centers_[:,0]);
   y = list(cluster.cluster_centers_[:,1]);
   z = list(cluster.cluster_centers_[:,2]);
   c = list(cluster.labels_);

   scipy.io.savemat('k_labels.mat', {'labels':cluster.labels_})
   scipy.io.savemat('k_coords.mat', {'coords':cluster.cluster_centers_})

   # plot the results of kmeans
   cmap = colors.Colormap('hot');
   norm  = colors.Normalize(vmin=1, vmax=10);
   s = 64;
   fig = plt.figure();
   ax = fig.add_subplot(111,projection='3d');
   Axes3D.scatter3D(ax,x,y,z,s=s);
   plt.show(fig);
   
   return cluster.cluster_centers_,cluster.labels_;