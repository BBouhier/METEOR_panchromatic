import numpy as np
from astropy.io import fits as pyfits
from astropy.table import Table, Column
from pathlib import Path
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
    p = Path(folder_path)
    fits = p.glob('*_anchored_star_galaxy_corrected.fits')
    txt = p.glob('*.txt')
    return list(fits), list(txt)

def create_catalogue(folder):
    fits,txt_ellips = get_files_in_folder(folder)
    ellips = Table.read(txt_ellips[0],format='ascii')
    a = ellips["a"]
    b = ellips["b"]
    y_source = ellips["y_centre"]
    x_source = ellips["x_centre"]
    theta = ellips["theta"]
    
    cube = np.array([pyfits.open(im)[0].data for im in fits])

    x,y = np.meshgrid(np.arange(np.shape(cube)[1]), np.arange(np.shape(cube)[2]))
    d = (np.cos(theta)*(x_source-x)+np.sin(theta)*(y_source-y))**2/a**2 + (np.sin(theta)*(x_source-x)-np.cos(theta)*(y_source-y))**2/b**2
    w = np.where((d < 1) & np.isfinite(np.max(cube, axis=0)) & (np.min(cube, axis=0) > 0))
    
    t = Table()
    t.add_column(Column([f"{c1}-{c2}" for c1,c2 in zip(x[w],y[w])],name = "id")) 
    t.add_column(Column(20.,name = "distance"))
    t.add_column(Column(0.,name = "redshift"))
    dic_lmbda_col = {'275': 'hst.wfc3.F275W',
                       "336":'hst.wfc3.F336W',
                       "438":'hst.wfc3.F438W',
                       "555": 'hst.wfc3.F555W',
                       "814": 'hst.wfc3.F814W',     
                       "200": 'jwst.nircam.F200W',
                       "300": 'jwst.nircam.F300M',
                       "335": 'jwst.nircam.F335M',
                       "360": 'jwst.nircam.F360M',
                       "770": 'jwst.miri.F770W',
                       "1000": 'jwst.miri.F1000W',
                       "1130": 'jwst.miri.F1130W',
                       "2100": 'jwst.miri.F2100W'}   
    
    for im_name in fits:
        lmbda = extract_lambda(str(im_name))
        image = pyfits.open(im_name)[0].data[w]
        t.add_column(Column(image,name = dic_lmbda_col[lmbda]))

    t.write('catalogg.txt', format='ascii', overwrite=True)
    t.write('cataloggue.fits', format='fits', overwrite=True)

if __name__ == '__main__':
    print('Hello there, here is a script to create a catalogue for JWST and HST fits files')
    print("Please enter the path to the directory containing the images in the different bands")
    folder_images = input("Path --> ")
    create_catalogue(folder_images)
    print("###############################################################################")
    print('The catalogue has been saved in the following folder: ')
    print(folder_images)
    print("###############################################################################")
