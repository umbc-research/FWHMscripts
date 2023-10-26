import astropy
import scipy
import numpy
import astropy
from astropy.io import fits
from pathlib import Path
from matplotlib import pyplot as plt
import profileFitting as pffit 


"""
GOALs:
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
TO DOs:
take userinput for which fits file
export images
test profiles compared to other known scripts
figure out how to display data, histrogram from pychart?
optimize
"""

fits_filename=None
hdul=None

#function to load in fits data, will eventually give a list of options and ask for num 1-9
#TODO:make auto var a boolean, if true then just find the most recently edited fits
def loadFits(auto=False):
  global fits_filename

  fits_filename = input("Name of .fits, exlcuding .fits:\t")
  if(fits_filename==""):
    fits_filename='theta_herc_00001'
  
  global hdul
  hdul = fits.open('{}{}{}{}'.format(Path.cwd(),'/fitsfiles/',fits_filename,'.fits'))


#function to create a profile fit for fits files, 0=vert, 1=horz, 2=2D, 3=radial
#TODO: litterally the entire thing, export images to somewhere? possibly make gui display for this whole thing?
def profileFit(type={0,1,2,3}):
  if(type==0 or type==1):

    return None
  elif(type==2):

    #recursively call for 0 and 1 and combine to get 2d graph

    return None
  else:
    
    return None


#function to calculate FWHM for a profile fit, returns float value(may change this) of the FWHM data 
#returns value in arcseconds, not pixels
#TODO:make this into a live enoivrment thing which can dynamically show a fwhm for different fits or something like that?
def fwhmCalc(prfFit):

  return None



#TODO: Design workflow for 'batch' FITS file FWHM identification.
#  - Read in (as command-line-input) the name of a directory and operate on all FITS files within
if __name__ == '__main__':
  print("work in progress as of 10/23, contact Olivia Chiarini (c241@umbc.edu) with questions")
  
  loadFits()
  
  done=False
  while(not done):
    try:
      typeFit=int(input('What type of fit? 0 for vert, 1 for horz, 2 for 2d, 3 for radial, 4 for all and compare, 5 for esc:\t'))
      match typeFit:
        case 0:
          pffit.fitCase0(hdul)
        case 1:
          pffit.fitCase1(hdul)
        case 2:
          pffit.fitCase2(hdul)
        case 3:
          pffit.fitCase3(hdul)
        case 4:
          pffit.fitCase4(hdul)  
        case _:
          exit()
      done=True
    except(Exception):
     print("Make sure your input was a number, and your fits file is correct!")

  #hdul.info()

  hdul.close()