import astropy
import scipy
import numpy
import astropy.io

"""
Authors:Olivia Chiarini
Research Advisors: Dr. Meyer and Roy Prouty
Purpose: To analyze FWHM data from .fits files given by UMBC's Obsertvatory, 
compare these to SHARPCAP to see how they are getting their histograms 
for FWHM and focusing analyizse. 
"""

"""
GOALS:
Load in FITS file
Find Sources
Fit Profiles
  -Vertical
  -Horizontal
  -Radial
 Extract STDV
Calculate FWHM
Convert from pixels at arcseconds
"""

"""
TO DOS:
take userinput
make folder for users to place fits files in
make filepath dynamio
figure out how to display data
"""

if __name__ == '__main__':
  print("work in progress sorry :(")
  



#STEP ONE: load in fits file
def loadFits():
  fits_filename = fits.util.get_testdata_filepath('test0.fits')
  hdul = fits.open(fits_filename)


#STEP TWO: find sources


#STEP THREE:fit profiles


#STEP FOUR: extract STDVS


#STEP FIVE: convert from pixels ot arcseconds


#STEP SIX: display and export this info


