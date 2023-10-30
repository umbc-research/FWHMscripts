import scipy
import numpy as np
  #Hi Olivia 
import astropy
from astropy.io import fits
from photutils.detection import DAOStarFinder,IRAFStarFinder

from pathlib import Path

from matplotlib import pyplot as plt
from sys import argv

#Author-defined pckgs
import profileFitting as pffit

"""
OVERALL TODOS:
  - Design workflow for 'batch' FITS file FWHM identification.
  - Read in (as command-line-input) the name of a directory and operate on all FITS files within
  - Some of the parameters in starfind should be made to update based on the measures of central tendancy above
  - Adapt source code to work for all sources within some brightness threshold
  - Build in support for changing plate scale (0.0317) for focal lengths
      config file?
  - Debug charts with curves displayed as well as original data
  - Find a  way to extract bounds from the FITS file or make reasonable guesses
    Roy suggests supplying the function with bounds based on FITS data/physical constraints of our detectors    
     
"""

fits_filename=None
hdul=None

#function to load in fits data
def loadFits(fits_filename):
  return fits.open(f'{fits_filename}')[0]

#function to calculate FWHM for a profile fit, returns float value(may change this) of the FWHM data 
#returns value in arcseconds, not pixels

if __name__ == '__main__':
  print("work in progress as of 10/23, contact Olivia Chiarini (c241@umbc.edu) with questions")

  try:
    fits_filename = argv[1]
  except IndexError:
    print("No file path provided as argument to python script.\nExiting")
    exit(0)

  ####Setting size of subFrames
  # leave as int!
  sF_length = 50
  
  hdul = loadFits(fits_filename)

  #Assumed to be in um
  pixSize = float(hdul.header['XPIXSZ'])

  data = hdul.data

  ####Find sources
  mean, median, std, max = np.mean(data), np.median(data), np.std(data), np.max(data)

  starFind = DAOStarFinder(threshold=median, fwhm=20.0, sky=500, exclude_border=True, brightest=10, peakmax=70000)
  sourcesList = starFind(data)
  
  ####Extract a subframe
  #Find center of first (brightest) source
  #  This code is not perfect, but works for a single source (here the 0th)
  xC, yC = sourcesList[0]['xcentroid'], sourcesList[0]['ycentroid']

  #Slicing (and matrices in general in python) are (row, col)
  subFrame = data[int(yC - sF_length):int(yC + sF_length),int(xC - sF_length):int(xC + sF_length)]
  
  ####Extract profile data (ideally centered on sources)
  
  ##Extract Radial Profile Data
  # Assuming subFrame is square
  radial_data_raw = pffit.extract_radial_data(subFrame, xC=sF_length, yC=sF_length)[:sF_length]
  radial_profile = np.concatenate((radial_data_raw[::-1], radial_data_raw))
  
  ##Extract Horizontal Data
  horiz_profile = subFrame[sF_length,:]

  ##Extract Vertical Data
  verti_data = subFrame[:,sF_length]

  ####Perform Fits
  x_data = np.linspace(0, sF_length*2, sF_length*2)

  radialParams = pffit.fit_gaussian_1d(x_data,radial_profile)

  horizParams = pffit.fit_gaussian_1d(x_data,horiz_profile)

  vertiParams = pffit.fit_gaussian_1d(x_data,verti_data)

  ####Convert the STD to FWHM and convert to " (arcsec) based on FITS header
  radFWHM = 2.355*radialParams[1] * 0.0317 * pixSize
  horizFWHM = 2.355*horizParams[1] * 0.0317 * pixSize
  vertiFWHM = 2.355*vertiParams[1] * 0.0317 * pixSize

  print(radFWHM, horizFWHM, vertiFWHM)

  ####Generate Plots, if desired

