#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:09:17 2023

@author: bbouhier

This code is to patch unwanted stars and galaxies on images of a galaxie. This code can treat HST and JWST fits files.
First step is to open each images with DS9 in order to write in a text file the different parameters of the ellipses
to patch galaxies. The stars will be patched automaticcaly with the gaia query


"""
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.table import Table
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import os 
import re

plt.close('all')

def gaia(filename):
    hst = re.search(r'sci', filename)
    jwst = re.search(r'anchored', filename)    
    if hst:
        RA = fits.open(filename)[0].header["RA_TARG"]
        DEC = fits.open(filename)[0].header["DEC_TARG"]
    elif jwst:
        RA = fits.open(filename)[1].header["RA_V1"]
        DEC = fits.open(filename)[1].header["DEC_V1"] 
    coord = SkyCoord(ra=RA, dec=DEC, unit=( u.degree, u.degree), frame='icrs', equinox = 'J2000')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(0.039, u.degree))
    r = j.get_results()
    RA_stars = r["ra"].data
    DEC_stars = r["dec"].data
    parallax_stars = r["parallax"].data
    pmra_stars = r["pmra"].data
    pmdec_stars = r["pmdec"].data
    return RA_stars, DEC_stars, parallax_stars, pmra_stars, pmdec_stars

def get_files_in_folder(folder_path):
    folder_path = os.path.abspath(folder_path)

    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        files_in_folder = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, file))]
        return files_in_folder
    else:
        print(f"Error: folder '{folder_path}' does not exist.")
        return None
    

def create_folder(folder):
    first_path = folder[0]
    galaxy = os.path.dirname(first_path)
    parent_directory_path = os.path.dirname(os.path.dirname(first_path))
    new_folder_name = os.path.basename(galaxy) + '_patch_correctionned'
    new_folder_path = os.path.join(parent_directory_path, new_folder_name)

    if os.path.exists(new_folder_path):
        print("Path", new_folder_path, 'already exists. Folder not created.')
        return False, new_folder_path
    
    os.makedirs(new_folder_path)
    print("New folder created:", new_folder_path) 
    return True, new_folder_path

def extract_lambda(filename):
    match = re.search(r'f(\d+)([a-zA-Z])', filename, re.IGNORECASE)
    if match:
        lambdaa = float(match.group(1))
        largeur = match.group(2)
        return lambdaa, largeur      
    else:
        print('No wavelength in title for {}'.format(filename))
        return None

def patch_stars(folder, plot_result = True):
    """
    Parameters
    ----------
    folder: path directory with all the jwst and hst image of the same galaxy, and the files for galaxies correction
    Returns
    -------
    None.

    """
    folder = get_files_in_folder(folder)
    
    created_folder, new_folder_path = create_folder(folder)
    
    txt_galaxies = next((file for file in folder if "txt_galaxies" in file), None)
    
    if not txt_galaxies:
        print("Error: missing .txt file to patch galaxies")
        
    data = Table.read(txt_galaxies, format='ascii')
    
    query = gaia(folder[0])
    
    ra_stars, dec_stars, parallax_stars, pmra_stars, pmdec_stars = query[0], query[1], query[2], query[3], query[4]

    ra_new = []
    dec_new = []
    threshold = 4.8e-8
    r = 10 #radius of the patch for stars
    for l in range (len(parallax_stars)):
        if parallax_stars[l]>threshold:
            ra_new.append(ra_stars[l])
            dec_new.append(dec_stars[l])
           
    for filename in folder:
        jwst = re.search(r'anchored', filename)
        hst = re.search(r'sci', filename)
        text = re.search(r'.txt', filename)
        output_filename = os.path.join(new_folder_path, os.path.basename(filename).replace('.fits', '_star_galaxy_corrected.fits'))
        
        if hst:
            hdu = fits.open(filename)[0]
            image_open = hdu.data
            hdul = hdu.header
            wavelength = extract_lambda(filename)
            lambdaa = wavelength[0]
            
            a_ell_arsec = []
            b_ell_arsec = []
            
            for a_ell in data['a']:
                convert_a = a_ell/(hdul['CD2_2']*3600)
                a_ell_arsec.append(convert_a)
                
            for b_ell in data['b']:
                convert_b = b_ell/(hdul['CD2_2']*3600)
                b_ell_arsec.append(convert_b)
            
            ra_g, dec_g, a_g, b_g, theta_g = data['ra'], data['dec'], a_ell_arsec, b_ell_arsec, data['theta']
            theta_rad = np.radians(theta_g)
            wcs = WCS(hdu.header)
            
            for l in range (len(ra_new)): 
                x_source, y_source = wcs.wcs_world2pix(ra_new[l], dec_new[l], 1)
                x_source= int(x_source)
                y_source = int(y_source)
                length_rec_half = int(3*r)
                width_rec_half = int(3*r)
                noise_array = []
            
                for i in range (x_source-length_rec_half,x_source+length_rec_half):
                    for j in range (y_source-width_rec_half,y_source+width_rec_half):
                        if (i-x_source)**2/r**2 + (j-y_source)**2/r**2 >1:
                            noise_array.append(image_open[i,j])

                mean = np.mean(noise_array)
                std = np.std(noise_array)
                
                for i in range (x_source-length_rec_half,x_source+length_rec_half):
                    for j in range (y_source-width_rec_half,y_source+width_rec_half):    
                        if (i-x_source)**2/r**2 + (j-y_source)**2/r**2 <1:
                            image_open[i,j]= np.random.normal(mean,std)

            for l in range (len(ra_g)): 
                y_source, x_source = wcs.wcs_world2pix(ra_g[l], dec_g[l], 1)
                x_source = int(x_source)
                y_source = int(y_source)
                length_rec_half = int(3*a_g[l])
                width_rec_half = int(3*b_g[l])
                noise_array = []
            
                for i in range (x_source-length_rec_half,x_source+length_rec_half):
                    for j in range (y_source-width_rec_half,y_source+width_rec_half):
                        if ((np.cos(theta_rad[l])*(x_source-i)+np.sin(theta_rad[l])*(y_source-j))**2/a_g[l]**2 
                            + (np.sin(theta_rad[l])*(x_source-i)-np.cos(theta_rad[l])*(y_source-j))**2/b_g[l]**2) >1:
                            noise_array.append(image_open[i,j])

                mean = np.mean(noise_array)
                std = np.std(noise_array)
                
                for i in range (x_source-length_rec_half,x_source+length_rec_half):
                    for j in range (y_source-width_rec_half,y_source+width_rec_half):
                        if ((np.cos(theta_rad[l])*(x_source-i)+np.sin(theta_rad[l])*(y_source-j))**2/a_g[l]**2 
                            + (np.sin(theta_rad[l])*(x_source-i)-np.cos(theta_rad[l])*(y_source-j))**2/b_g[l]**2) <1:
                            image_open[i,j]= np.random.normal(mean,std)
        
            fits.writeto(output_filename, image_open, header=fits.getheader(filename), overwrite=True)
            print('New file created:', output_filename)
            
        elif jwst:
            hdu = fits.open(filename)[1]
            image_open = hdu.data
            hdul = hdu.header
            wavelength = extract_lambda(filename)
            lambdaa = wavelength[0]
            
            a_ell_arsec = []
            b_ell_arsec = []
            
            for a_ell in data['a']:
                convert_a = a_ell/(hdul['CDELT1']*3600)
                a_ell_arsec.append(convert_a)
                
            for b_ell in data['b']:
                convert_b = b_ell/(hdul['CDELT1']*3600)
                b_ell_arsec.append(convert_b)
            
            ra_g, dec_g, a_g, b_g, theta_g = data['ra'], data['dec'], a_ell_arsec, b_ell_arsec, data['theta']
            theta_rad = np.radians(theta_g)
            wcs = WCS(hdu.header)

            for l in range(len(ra_g)):
                y_source, x_source = wcs.wcs_world2pix(ra_g[l], dec_g[l], 1)
                x_source = int(x_source)
                y_source = int(y_source)
                length_rec_half = int(3 * a_g[l])
                width_rec_half = int(3 * b_g[l])
                noise_array = []
            
                pixel_value = image_open[y_source, x_source]
            
                if not np.isnan(pixel_value):
                    for i in range(x_source - length_rec_half, x_source + length_rec_half):
                        for j in range(y_source - width_rec_half, y_source + width_rec_half):
                            if ((np.cos(theta_rad[l]) * (x_source - i) + np.sin(theta_rad[l]) * (y_source - j))**2 / a_g[l]**2 +
                                    (np.sin(theta_rad[l]) * (x_source - i) - np.cos(theta_rad[l]) * (y_source - j))**2 / b_g[l]**2) > 1:
                                pixel = image_open[i, j]
                                if not np.isnan(pixel):
                                    noise_array.append(pixel)
            
                    mean = np.mean(noise_array)
                    std = np.std(noise_array)
            
                    for i in range(x_source - length_rec_half, x_source + length_rec_half):
                        for j in range(y_source - width_rec_half, y_source + width_rec_half):
                            if ((np.cos(theta_rad[l]) * (x_source - i) + np.sin(theta_rad[l]) * (y_source - j))**2 / a_g[l]**2 +
                                    (np.sin(theta_rad[l]) * (x_source - i) - np.cos(theta_rad[l]) * (y_source - j))**2 / b_g[l]**2) < 1 and not np.isnan(image_open[i, j]):
                                image_open[i, j] = np.random.normal(mean, std)
            
            fits.writeto(output_filename, image_open, header=fits.getheader(filename), overwrite=True)
            print('New file created:', output_filename)
            
        elif text:
            pass
        
        else:
            print(f"Error: Unknown instrument type for file {filename}. It should be either HST or JWST.")
                
        if plot_result:
            if jwst:
                fig, ax = plt.subplots(subplot_kw={'projection': wcs}, figsize=(8, 8))
                ax.imshow(image_open)
            
                ra_new = []
                dec_new = []
                for l in range(len(parallax_stars)):
                    if parallax_stars[l] > threshold:
                        ra_new.append(ra_stars[l])
                        dec_new.append(dec_stars[l])
            
                x, y = wcs.wcs_world2pix(ra_g, dec_g, 1)
                ax.scatter(x, y, marker='D', color='Green', label='Galaxies')
                ax.legend()
                plt.title(f'JWST at {lambdaa} nm')
                plot_filename = os.path.join(new_folder_path, f'JWST_{lambdaa}_plot.png')
                plt.savefig(plot_filename)  
                plt.close('all')
            
            elif hst:
                vmin = np.percentile(image_open, 5)
                vmax = np.percentile(image_open, 99.95)
                
                fig, ax = plt.subplots(subplot_kw={'projection': wcs}, figsize=(8, 8))
                ax.imshow(image_open, vmin=vmin, vmax=vmax)
            
                x, y = wcs.wcs_world2pix(ra_new, dec_new, 1)
                ax.scatter(y, x, marker='+', color='red', label='Stars')
                x, y = wcs.wcs_world2pix(ra_g, dec_g, 1)
                ax.scatter(x, y, marker='D', color='Green', label='Galaxies')
                ax.legend()
                plt.title(f'HST at {lambdaa} nm')   
                plot_filename = os.path.join(new_folder_path, f'HST_{lambdaa}_plot.png')
                plt.savefig(plot_filename) 
                plt.close('all')
                
    print("###############################################################################")
    print('All the corrected fits files and plots have been save in the following folder: ')
    print(new_folder_path)
    print("###############################################################################")


if __name__ == "__main__":
    print('Hello there, here is a script to correct JWST and HST fits files')
    print("Please enter the path to the directory with all the files for one galaxy")
    folder = input("Path --> ")
    patch_stars(folder)
