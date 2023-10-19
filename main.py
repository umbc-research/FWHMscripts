import astropy
import scipy
import numpy
import astropy
from astropy.io import fits
from pathlib import Path

"""
Authors:Olivia Chiarini
Research Advisors: Dr. Meyer and Roy Prouty
Purpose: To analyze FWHM data from .fits files given by UMBC's Obsertvatory, 
compare these to SHARPCAP to see how they are getting their histograms 
for FWHM and focusing analyizse. 
"""

"""
GOALS:
Load in FITS file [x]
Find Sources[ ]
Fit Profiles[ ]
  -Vertical[ ]
  -Horizontal[ ]
  -Radial[ ]
 Extract STDV[ ]
Calculate FWHM[ ]
Convert from pixels at arcseconds[ ]
"""

"""
TO DOS:
take userinput for which fits file
export images
test profiles compared to other known scripts
figure out how to display data, histrogram from pychart?
optimize
"""



fits_filename=None
hdul=None


def loadFits(fitsname):
  global hdul
  hdul = fits.open('{}{}{}{}'.format(Path.cwd(),'/fitsfiles/',fitsname,'.fits'))


if __name__ == '__main__':
  print("work in progress as of 10/23, contact Olivia Chiarini (c241@umbc.edu) with questions")
  fits_filename = input("Name of .fits, exlcuding .fits:\t")
  if(fits_filename==""):
    fits_filename='theta_herc_00001'

  loadFits(fits_filename)
  hdul.info()

  hdul.close()