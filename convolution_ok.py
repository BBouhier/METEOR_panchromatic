# -*- coding: utf-8 -*-
# Author: Héctor Salas

import numpy as np
from astropy.io import fits
from scipy.ndimage import zoom
from astropy.convolution import convolve_fft
import astropy
import readline
readline.parse_and_bind('tab: complete')
readline.set_completer_delims(' \t\n')
from astropy.wcs import WCS
import os 
import re

#complete this function to get a more invormative header
def update_header(header_i, header_k):
    header = header_i.copy()
    header['history'] = ('Convolved version of created with convolve_images.py')
    return header

def do_the_convolution(image, image_h, kernel, kernel_h):
    """Function that convolve an image with a kernel
    Inputs:
        image:  image data
        image_h: image header
        kernel: kernel data
        kernel_h: kernel_header
    Outputs:
    """

    # Get image and kernel pixel_scale
    wcs_image_h = WCS(image_h)
    wcs_kernel_h = WCS(kernel_h)

    pixel_scale_i = astropy.wcs.utils.proj_plane_pixel_scales(wcs_image_h)[0]
    pixel_scale_k = astropy.wcs.utils.proj_plane_pixel_scales(wcs_kernel_h)[0]
    # resize kernel if necessary
    if pixel_scale_k != pixel_scale_i:
        ratio = pixel_scale_k / pixel_scale_i
        size = ratio*kernel.shape[0]
        # ensure a odd kernel
        if round(size) % 2 == 0:
            size += 1
            ratio = size / kernel.shape[0]
        kernel = zoom(kernel, ratio) / ratio**2
        
    # do convolution
    if len(np.shape(image)) == 2:
        convolved = convolve_fft(image, kernel, nan_treatment='interpolate',
                             normalize_kernel=False, preserve_nan=True,
                             boundary='fill', fill_value=0., allow_huge=True)
    else:
        convolved = []
        for i in range(len(np.shape(image))+1):
            convolved_i = convolve_fft(image[i], kernel,
                                       nan_treatment='interpolate',
                                       normalize_kernel=False,
                                       preserve_nan=True, fft_pad=True,
                                       boundary='fill', fill_value=0., allow_huge=True)
            convolved.append(list(convolved_i))
        convolved = np.asarray(convolved)

    return convolved, kernel

def extract_lambda(filename):
    match = re.search(r'f(\d+)([a-zA-Z])', filename, re.IGNORECASE)
    if match:
        lambdaa = str(match.group(1))
        return lambdaa     
    else:
        print('No wavelength in title for {}'.format(filename))
        return None
    
def extract_size_k(filename):
    match = re.search(r's(\d+\.\d+)', filename, re.IGNORECASE)
    match2 = re.search(r's(\d+)', filename, re.IGNORECASE)
    if match:
        size = str(match.group(1))
        return size
    elif match2:
        size = str(match2.group(1))
        return size
    else:
        print('No wavelength in title for {}'.format(filename))
        return None
    
def get_files_in_folder(folder_path):
    folder_path = os.path.abspath(folder_path)

    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        files_in_folder = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, file))]
        return files_in_folder
    else:
        print(f"Error: folder '{folder_path}' does not exist.")
        return None
    
def get_fits_files_in_folder(folder_path):
    folder_path = os.path.abspath(folder_path)

    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        fits_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, file)) and file.lower().endswith('.fits')]
        return fits_files
    else:
        print(f"Error: folder '{folder_path}' does not exist.")
        return None

def create_folder(folder,gaussian_size):
    first_path = folder[0]
    galaxy = os.path.dirname(first_path)
    parent_directory_path = os.path.dirname(os.path.dirname(first_path))
    new_folder_name = os.path.basename(galaxy) + '_convolved_g{}'.format(gaussian_size)
    new_folder_path = os.path.join(parent_directory_path, new_folder_name)

    if os.path.exists(new_folder_path):
        print("Path", new_folder_path, 'already exists. Folder not created.')
        return False, new_folder_path
    
    os.makedirs(new_folder_path)
    print("New folder created:", new_folder_path) 
    return True, new_folder_path

if __name__ == '__main__':
    print('Hello there, here is a script to convolve JWST and HST fits files')
    print("Please enter the path to the directory with all the files to convolve")
    folder_images = input("Path --> ")
    print("Please enter the path to the directory with all the kernels")
    folder_kernels = input("Path --> ")
    print("Please choose a Gaussian size between: 0.85, 1.15, 4, 7.5, 15")
    gaussian_size = input("Gaussian size --> ")

    if gaussian_size not in ['0.85', '1.15', '4', '7.5', '15']:
        print("invalid Gaussian size ")
    
    folder_kernels = get_files_in_folder(folder_kernels)
    folder_images = get_fits_files_in_folder(folder_images)
    created_folder, new_folder_path = create_folder(folder_images,gaussian_size)
    
    selected_kernels=[]
    
    #extract size from kernels
    for kernel in folder_kernels:
        size = extract_size_k(kernel)
        if size == gaussian_size:
            selected_kernels.append(kernel)
    
    for image in folder_images:
        hst = re.search(r'sci', image)
        jwst = re.search(r'anchored', image)
        lmbda_im = extract_lambda(image)
        if jwst:
            for k in selected_kernels:
                wl = extract_lambda(k)
                if wl == lmbda_im:
                    kernel_name = k
                    print(kernel_name)
                    break

            kernel =  fits.open(kernel_name)[0].data
            header_k = fits.open(kernel_name)[0].header
            image_data =  fits.open(image)[0].data
            header_i = fits.open(image)[0].header

        elif hst:
            for k in selected_kernels:
                wl = extract_lambda(k)
                if wl == '200':
                    kernel_name = k
                    break
            
            kernel =  fits.open(kernel_name)[0].data
            header_k = fits.open(kernel_name)[0].header
            image_data =  fits.open(image)[0].data
            header_i = fits.open(image)[0].header
        
        kernel = kernel/np.max(kernel)

        convolved, kernel = do_the_convolution(image_data, header_i, kernel, header_k)
        
        header_new = update_header(header_i, header_k)
        convolved_name = os.path.join(new_folder_path, f"convolved_g{gaussian_size}_{os.path.basename(image)}")
        fits.writeto(convolved_name, convolved, header_new, overwrite=True)
    
    print("###############################################################################")
    print('All the corrected fits files and plots have been save in the following folder: ')
    print(new_folder_path)
    print("###############################################################################")