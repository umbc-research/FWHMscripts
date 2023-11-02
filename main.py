import scipy
import numpy as np

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
  - Work on dumping all findings and necessary info for reproducing to a log file. Save figure, too.  
"""

fits_filename=None
hdul=None

#function to load in fits data
def loadFits(fits_filename):
  return fits.open(f'{fits_filename}')[0]

#function to calculate FWHM for a profile fit, returns float value(may change this) of the FWHM data 
#returns value in arcseconds, not pixels



#takes in file path, and string, checks if string is debug
if __name__ == '__main__':
  print("work in progress as of 10/23, contact Olivia Chiarini (c241@umbc.edu) with questions")

  try:
    fits_filename = argv[1]
  except IndexError:
    print("No file path provided as argument to python script.\nExiting")
    exit(0)

  try:
    debugCheck=(argv[2]=="debug")
    print("Running in debug mode")
  except:
    pass

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
  radial_data = np.concatenate((radial_data_raw[::-1], radial_data_raw))
  
  ##Extract Horizontal Data
  horiz_data = subFrame[sF_length,:]

  ##Extract Vertical Data
  verti_data = subFrame[:,sF_length]

  ####Perform Fits
  x_data = np.linspace(0, sF_length*2, sF_length*2)

  radialParams = pffit.fit_gaussian_1d(x_data,radial_data)
  horizParams = pffit.fit_gaussian_1d(x_data,horiz_data)
  vertiParams = pffit.fit_gaussian_1d(x_data,verti_data)

  #####Get Fitted Profiles
  radial_fit = pffit.gaussian_1d(x_data, *radialParams)
  horiz_fit  = pffit.gaussian_1d(x_data, *horizParams)
  verti_fit  = pffit.gaussian_1d(x_data, *vertiParams)

  ##Calculate Residuals
  ## SQRT( SUM( SQUARED_DIFFERENCES  ) ) / numDataPts
  radialResidual = np.sqrt( sum( (radial_fit - radial_data) **2 ) ) / (2*sF_length)
  horizResidual  = np.sqrt( sum( (horiz_fit  - horiz_data)  **2 ) ) / (2*sF_length)
  vertiResidual  = np.sqrt( sum( (verti_fit  - verti_data)  **2 ) ) / (2*sF_length)

  print(f"Fits completed with the following residuals\nRadial: {radialResidual:0.3f}\nHorizontal: {horizResidual:0.3f}\nVertical: {vertiResidual:0.3f}\n")

  ####Convert the STD to FWHM and convert to " (arcsec) based on FITS header
  radFWHM   = 2.355*radialParams[1] * 0.0317 * pixSize
  horizFWHM = 2.355*horizParams[1]  * 0.0317 * pixSize
  vertiFWHM = 2.355*vertiParams[1]  * 0.0317 * pixSize

  print(f"Radial FWHM: {radFWHM:0.3f}\nHorizontal FWHM: {horizFWHM:0.3f}\nVertical FWHM: {vertiFWHM:0.3f}")
  
  ####Generate Plots

  fig,charts =plt.subplots(2,2, figsize=(10,8))
  
  charts[0,0].plot(x_data, horiz_data, 'ko', markersize=2)
  charts[0,0].plot(x_data, horiz_fit,'tab:blue', linestyle='dashed',)
  charts[0,0].set_title("Horizontal Fit")
  charts[0,0].set(xlabel='Pixel in SubFrame',ylabel='Counts')
  charts[0,0].grid(1)

  charts[0,1].plot(x_data, verti_data, 'ko', markersize=2)
  charts[0,1].plot(x_data, verti_fit,'tab:red', linestyle='dashed',)
  charts[0,1].set_title("Vertical Fit")
  charts[0,1].set(xlabel='Pixel in SubFrame',ylabel='Counts')
  charts[0,1].grid(1)

  charts[1,0].plot(x_data,radial_data, 'ko', markersize=2)
  charts[1,0].plot(x_data,radial_fit,'tab:purple', linestyle='dashed')
  charts[1,0].set_title("Radial Fit")
  charts[1,0].set(xlabel='Pixel in SubFrame',ylabel='Counts')
  charts[1,0].grid(1)

  charts[1,1].imshow(subFrame, cmap='gray')
  charts[1,1].plot(sF_length,sF_length, 'rx')
  charts[1,1].set_title("Source SubFrame")

  plt.suptitle("FWHM Curve Fitting")
  plt.tight_layout()
  plt.show()