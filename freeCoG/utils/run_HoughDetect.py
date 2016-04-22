import cv2
import cv2.cv as cv
import scipy
import nibabel as nib
import numpy as np
import os


def CT_hough(imgs_dir,img):

    # nav to dir with images to process
    os.chdir(imgs_dir);

    # load in thresholded and smoothed CT and extract data matrix
    im_l = nib.load(img); 
    img_data = im_l.get_data();

    # convert CT data to uint8 for opencv proc
    #cv_img = img_data.astype('uint8');

    # get iterable list of slices
    z_len = img_data[1,1,:].shape[0];# of slices in the z-direction
    x = np.linspace(1,z_len, z_len);
    x.astype('int');
    x = x[1:(len(x)-1)];


    # iterate through slices and apply hough transform to detect circles
    elect_coords = np.array([0, 0, 0]);
    for i in x:

        # slice through the z dimension of the CT
        slice = img_data[:,:,int(i)];

        #save the slice as .png
        scipy.misc.imsave('CT_slice.png', slice);

        #read in the save .png im slice
        ct_img = cv2.imread('CT_slice.png',0);

        # Convert image slice to grayscale for viewing
        cimg = cv2.cvtColor(ct_img,cv2.COLOR_GRAY2BGR);

        # run open-cv hough cirlces detection using gradient method
        circles = cv2.HoughCircles(ct_img,cv.CV_HOUGH_GRADIENT,1,4, param1=7,param2=5,minRadius=1,maxRadius=4);
    
    
        if (circles is None) == True:
            None;
        else:
            circles = circles[0,:,:];
            circles[:,2] = i;
            stackem = (elect_coords, circles);
            elect_coords = np.vstack(stackem);
        
    # throw out first row of zeros
    elect_coords = elect_coords[1:elect_coords[:,1].size,:];
    return elect_coords;
