# img_Convolve.py
# performs convoltion of a 3d img kernel over 
# an input target 3d image using the fourier
# image convolution method: inverse_fft(fft(img) .* fft(kernel))

import os
import scipy
import numpy
import nibabel

def img_conv(img,kernel):
    
    # use fsl to obtain fft of input img
    new_fname = 'fft' + img;
    os.system('fslfft %s %s' %(img, new_imgfname);
    
    # use fsl to obtain fft of input img
    new_fname = 'fft' + kernel;
    os.system('fslfft %s %s' %(kernel, new_kernelfname);
    
    # load in complex fft results
    fft_img = nibabel.load(new_imgfname);
    fft_img_data = fft_img.get_data();
    
    fft_kernel = nibabel.load(new_kernelfname);
    fft_kernel_data = fft_kernel.get_data();
    
    # perform point wise multiplication of complex fft outputs
    fft_mulRes = fft_img_data *. fft_kernel_data;
    
    # save the smoothed thresholded ct
    hdr = fft_img.get_header();
    affine = fft_img.get_affine();
    N = nib.Nifti1Image(fft_mulRes,affine,hdr);
    new_fname = 'fft_mulRes.nii';
    N.to_filename(new_fname);
    
    # load kernel sphere
    kernel = nibabel.load('kernel_sphere.nii')
    kernel= kernel.get_data()

    # perform 3d image convolution using ndimage.convolve
    convolved = ndimage.convolve(CT_img_dat, kernel, mode ='constant', cval=0.0)
    
    #save convolved img to nifti
    N = nibabel.Nifti1Image(convolved,affine,hdr);
    N.to_filename(new_fname);


    
    
    
    