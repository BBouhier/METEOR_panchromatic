# METEOR_panchromatic
Fits_treatment

Hello there,

This code is a component of the METEOR project, where the goal is to develop a pipeline for correcting HST or JWST images of galaxies. The ultimate aim is to derive various characteristics of the galaxies.

As the project progresses, this readme file will be updated and expanded.

The first program, patch_stars_galaxies_HST_JWST.py, will correct automatically stars on JWST and HST images with a Gaia query. It can also correct  galaxies for HST files by using DS9 to select the best ellips to patch the galaxies.
You will need to write the Ra and Dec in degree of the center of the ellipses in a text file, with the a and b of the ellipses in arcseconds, and the tilt theta in degrees if needed.
So at the end you will have a txt file ra, dec, a, b, theta. This txt file needs to be in the same folder as the fits images of the galaxy you want to correct.

You will then give the folder with your images and this txt file as input to the program !
At the end, a new folder will be created with all the corrected images and the plots.

Secondly, the program convolution.py will convolve each corrected image with the same Kernel in order to have the same resolution for all. You will need to put your file of corrected images, your file of kernels and the size of the selected kernel as an input.
A new folder containing the convolved images will be created.

Then, with the program reprojection.py,  a new header, select will be created with the number of pixels you want, their size and position of the central pixel.
Thanks to that, you will be able to reproject your corrected and convolved images in the same grid, in order to have the same image size etc...

The last part of the correction process is to de-redden the images from the effect of the dust of the Milky Way. Again, you just have to give in input the path to the folder created by reprojection.py, and it will create a new folder. Note: only the HST images are de-redddened.

################
Catalog and CIGALE part
