#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 17:36:49 2023

@author: bbouhier
"""
import astropy
import reproject
import numpy as np
from astropy.io import fits
import os 
from pathlib import Path

def get_files_in_folder(folder_path):
    p = Path(folder_path)
    HST = p.glob('*sci_star_galaxy_corrected.fits')
    JWST = p.glob('*anchored_star_galaxy_corrected.fits')
    return list(HST), list(JWST)
    
def create_folder(folder):
    first_path = folder[0]
    galaxy = os.path.dirname(first_path)
    parent_directory_path = os.path.dirname(os.path.dirname(first_path))
    new_folder_name = os.path.basename(galaxy) + '_reprojected'
    new_folder_path = os.path.join(parent_directory_path, new_folder_name)

    if os.path.exists(new_folder_path):
        print("Path", new_folder_path, 'already exists. Folder not created.')
        return False, new_folder_path
    
    os.makedirs(new_folder_path)
    print("New folder created:", new_folder_path) 
    return True, new_folder_path

def create_footprint(folder):
    header = astropy.io.fits.Header(
        cards=["NAXIS", "NAXIS1", "NAXIS2", "CD1_1", "CD1_2", "CD2_1", "CD2_2", 
               "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",],copy=False)

    header["NAXIS"] = 2
    header["NAXIS1"] = header["NAXIS2"] = 1001
    header["CD1_1"] = header["CD2_2"] = (500/20E6) * 180 / np.pi  #degree/pix
    header["CD1_2"] = header["CD2_1"] = 0.0
    header["CRPIX1"] = 501
    header["CRPIX2"] = 501
    header["CRVAL1"] = 41.75039357231478 
    header["CRVAL2"] = -0.5694873788265216  
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    
    HST, JWST = get_files_in_folder(folder)
    created_folder, new_folder_path = create_folder(HST)
    
    for filename in HST:
        filename = str(filename)
        reprojected_image_name = os.path.join(new_folder_path, f"reprojected_{os.path.basename(filename)}")
        image = fits.open(filename)[0]  # HST
        reprojected_image, footprint = reproject.reproject_exact(image, header, reprojected_image_name)
        fits.writeto(reprojected_image_name, reprojected_image, header, overwrite=True)
            
    for filename in JWST:
        filename = str(filename)
        reprojected_image_name = os.path.join(new_folder_path, f"reprojected_{os.path.basename(filename)}")
        image = fits.open(filename)[0]  # JWST
        reprojected_image, footprint = reproject.reproject_exact(image, header, reprojected_image_name)
        fits.writeto(reprojected_image_name, reprojected_image, header, overwrite=True)

    print("###############################################################################")
    print('All the corrected fits files and plots have been save in the following folder: ')
    print(new_folder_path)
    print("###############################################################################")

if __name__ == "__main__":
    print('Hello there, here is a script to correct JWST and HST fits files')
    print("Please enter the path to the directory with all the files for one galaxy that are already corrected")
    folder = input("Path --> ")
    create_footprint(folder)