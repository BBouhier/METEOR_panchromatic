from dust_extinction.parameter_averages import O94
import numpy as np
from astroquery.ipac.irsa.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import os 
import re

def extract_lambda(filename):
    match = re.search(r'f(\d+)([a-zA-Z])', filename, re.IGNORECASE)
    if match:
        lambdaa = str(match.group(1))
        return lambdaa     
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
    
def create_folder(folder):
    first_path = folder[0]
    galaxy = os.path.dirname(first_path)
    parent_directory_path = os.path.dirname(os.path.dirname(first_path))
    new_folder_name = os.path.basename(galaxy) + '_dereddened'
    new_folder_path = os.path.join(parent_directory_path, new_folder_name)

    if os.path.exists(new_folder_path):
        print("Path", new_folder_path, 'already exists. Folder not created.')
        return False, new_folder_path
    
    os.makedirs(new_folder_path)
    print("New folder created:", new_folder_path) 
    return True, new_folder_path


def dereddening (image_name):

    dic_lmbda_pivot = {'275': 2709.7*u.AA,
                       "336": 3354.5*u.AA,
                       "438": 4326.2*u.AA,
                       "555": 5308.4*u.AA,
                       "814": 8039.1*u.AA}
    hst = re.search(r'sci', image_name)
    image = fits.open(image_name)[0].data
    if hst:
        hdu = fits.open(image_name)[0].header
        
        ra = hdu["CRVAL1"]         
        dec = hdu["CRVAL2"]
        coo = coord.SkyCoord(ra,dec, unit =u.deg, frame='fk4')#1385

        lmbda_im = extract_lambda(image_name)
        lmbda_pivot = dic_lmbda_pivot[lmbda_im]
        table = IrsaDust.get_query_table(coo,section='ebv')
        ebv = table['ext SFD min'][0]
        o94 = O94(Rv=3.1)
        ext = o94.extinguish(lmbda_pivot,Ebv=ebv)
        for i in range(len(image)):
            for j in range (len(image)):
                image[i,j]=image[i,j]/ext
                
    return image
        

if __name__ == '__main__':

    
    print('Hello there, here is a script to deredden HST fits files')
    print("Please enter the path to the directory with all the files to deredden")
    folder_images = input("Path --> ")


    folder_images = get_files_in_folder(folder_images)
    created_folder, new_folder_path = create_folder(folder_images)
        
    for image_name in folder_images:
        dereddened =  dereddening(image_name)
        header = fits.open(image_name)[0].header
        dereddened_name = os.path.join(new_folder_path, f"dereddened_{os.path.basename(image_name)}")
        fits.writeto(dereddened_name, dereddened, header, overwrite=True)
    
    print("###############################################################################")
    print('All the corrected fits files and plots have been save in the following folder: ')
    print(new_folder_path)
    print("###############################################################################")