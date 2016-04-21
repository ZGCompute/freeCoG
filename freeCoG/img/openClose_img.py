# openClose_img.py
# function performs binary opening and closing of an img 
# to get rid of small holes and white spots on a binary img


from scipy import ndimage
import nibabel as nib


def openClose_img(fname):
    
    # load in thresholded CT volume
    img = nib.load(fname);
    im = img.get_data();

    # perform binary open and closing to clean img
    open_img = ndimage.binary_opening(im);
    close_img = ndimage.binary_closing(open_img);

    # save binary opened/closed image
    hdr = img.get_header()
    affine = img.get_affine()
    N = nib.Nifti1Image(close_img,affine,hdr);
    new_fname = "closed_" + fname;
    N.to_filename(new_fname);
