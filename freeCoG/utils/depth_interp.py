# depth_interp.py
# function to interpolote depth electrode coordinates
# from two endpoints in fs surface RAS space

# Author: Zachary Greenberg
# Department of Neurological Surgery
# University of California, San Francisco
# Date Last Edited: January 26, 2015


import numpy as np
import scipy.io
import os

def interp_depth(elec1, elec2, numElecs, fname):

   # interpolate through x y and z from end to end
   xs = np.linspace(elec1[0], elec2[0], numElecs);
   ys = np.linspace(elec1[1], elec2[1], numElecs);
   zs = np.linspace(elec1[2], elec2[2], numElecs);

   # flip interp xyz coords to stack as one matrix
   xs = np.transpose(xs);
   ys = np.transpose(ys);
   zs = np.transpose(zs);
   stackem = (xs,ys,zs);
   elecmatrix = np.column_stack(stackem);

   #  save result
   scipy.io.savemat(fname, {'elecmatrix':elecmatrix});

