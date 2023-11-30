# METEOR_panchromatic
Fits_treatment

Hello,

This code is a component of the METEOR project, where the goal is to develop a pipeline for correcting HST or JWST images of galaxies. The ultimate aim is to derive various characteristics of the galaxies.

As the project progresses, this readme file will be updated and expanded.

The first program, patching.py, will correct automatically stars on JWST and HST images with a Gaia query. It can also correct 
galaxies for HST files by using DS9 to select the best ellipse to patch the galaxies.
You will need to write the Ra and Dec in degree of the center of the ellipse in a text file, with the a and b of the ellipses in arcsseconds, and the tilt theta in degree if needed.
So at th end you will have a txt file ra, dec, a, b, theta

Then the txt file needs to be in the same folder as the fits images of the galaxy you want to correct.

You will then give the folder with your images and this txt file as input to the program !
