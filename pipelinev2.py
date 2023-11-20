#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 15:24:30 2023

@author: bbouhier
"""

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

from astroquery.simbad import Simbad
from astroquery.ipac.ned import Ned
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
import glob
import os 
import astropy.io.fits as pyfits

plt.close('all')

jwst = ['/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F200W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F300W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F335W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F360W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F770W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F1000W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F1130W_anchored.fits',
        '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_F2100W_anchored.fits']

hst = ['/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_uvis_f275w_exp_drc_sci.fits',
       '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_uvis_f336w_exp_drc_sci.fits',
       '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_uvis_f438w_exp_drc_sci.fits',
       '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_uvis_f555w_exp_drc_sci.fits',
       '/home/bbouhier/Desktop/M2_S1/METEOR_PANCHROMATIC_EMISSION/data_one_drive_ngc_1087/NGC_1087/ngc1087_uvis_f814w_exp_drc_sci.fits',]

RA = fits.open(hst[0])[0].header["RA_TARG"]
DEC = fits.open(hst[0])[0].header["DEC_TARG"]
 
    
def gaia(ra,dec):
    coord = SkyCoord(ra=RA, dec=DEC, unit=( u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(0.039, u.degree))
    r = j.get_results()

    RA_stars = r["ra"].data
    DEC_stars = r["dec"].data
    parallax_stars = r["parallax"].data
    pmra_stars = r["pmra"].data
    pmdec_stars = r["pmdec"].data
    return RA_stars, DEC_stars, parallax_stars, pmra_stars, pmdec_stars

search = gaia(RA, DEC)


'''
Patch for stars
'''

def patch_stars (image, ra, dec, parallax, a, b, theta, pmra, pmdec, plot_result = False):#stars
    """
    Parameters
    ----------
    image : numpy.ndarray
        Image you want to modify (example: fits.open(file_name)[0].data)
    ra : np.array
        Array of the Right Ascension of the sources
    dec : np.array
        Array of the Declination of the sources
    a : int
        semi major axis of patch ellips
    b : int
        semi minor axis of patch ellips
    theta : in degree
        inclination of the ellips
    plot_result: false by default
        to plot the image corrected

    Returns
    -------
    None.

    """
    try:
        output_filename = os.path.join(os.path.dirname(image), os.path.basename(image).replace('.fits', '_star_corrected.fits'))
        hdu = fits.open(hst[0])[0]
        image_open = hdu.data
        wcs = WCS(hdu.header)
        theta = np.radians(theta)
        ra_new = []
        dec_new = []
        threshold = 4.8e-8
        for l in range (len(parallax)):
            if parallax[l]>threshold:
                ra_new.append(ra[l])
                dec_new.append(dec[l])
            # if pmra[l]>threshold or pmdec[l]>threshold:
            #     print(pmra[l],pmdec[l])
            #     ra_new.append(ra[l])
            #     dec_new.append(dec[l])
                
            
        for l in range (len(ra_new)):
            x_source, y_source = wcs.wcs_world2pix(ra_new[l], dec_new[l], 1)
            x_source= int(x_source)
            y_source = int(y_source)
            length_rec_half = int(3*a)
            width_rec_half = int(3*b)
            noise_array = []
        
            for i in range (x_source-length_rec_half,x_source+length_rec_half):
                for j in range (y_source-width_rec_half,y_source+width_rec_half):
                    if (np.cos(theta)*(i-x_source)-np.sin(theta)*(j-y_source))**2/a**2 + (np.sin(theta)*(i-x_source)+np.cos(theta)*(j-y_source))**2/b**2 >1:
                        noise_array.append(image_open[i,j])
         
            mean = np.mean(noise_array)
            std = np.std(noise_array)
            
            for i in range (x_source-length_rec_half,x_source+length_rec_half):
                for j in range (y_source-width_rec_half,y_source+width_rec_half):
                    if (np.cos(theta)*(i-x_source)-np.sin(theta)*(j-y_source))**2/a**2 + (np.sin(theta)*(i-x_source)+np.cos(theta)*(j-y_source))**2/b**2 <1:
                        image_open[i,j]= np.random.normal(mean,std)
                        
        fits.writeto(output_filename, image_open, header = fits.getheader(image), overwrite=True)
        print('New file created:', output_filename)
    except Exception as e:
        print('Error:', str(e))
        
        if plot_result:
            
            vmin = np.percentile(image_open, 5)
            vmax = np.percentile(image_open, 99.95)
        
            fig, ax = plt.subplots(subplot_kw={'projection': wcs}, figsize=(8, 8))
            im = ax.imshow(image_open, vmin=vmin, vmax=vmax)
        
            ra_new = []
            dec_new = []
            for l in range(len(search[2])):
                if search[2][l] > threshold:
                    ra_new.append(search[0][l])
                    dec_new.append(search[1][l])
        
            ax.coords.grid(True, color='white', ls='solid')
            ax.coords[0].set_axislabel('Galactic Longitude')
            ax.coords[1].set_axislabel('Galactic Latitude')
        
            overlay = ax.get_coords_overlay('fk5')
            overlay.grid(color='black', ls='dotted')
            overlay[0].set_axislabel('Right Ascension (J2000)')
            overlay[1].set_axislabel('Declination (J2000)')
        
            x, y = wcs.wcs_world2pix(ra_new, dec_new, 1)
            ax.scatter(y, x, marker='o', color='red', label='Sources')
            ax.legend()
            plt.show()


a_stars = 15
b_stars = a_stars

for i, filename in enumerate(hst):
    star_corr = patch_stars(filename, search[0], search[1], search[2], a_stars, b_stars, 30, search[3], search[4])
    
