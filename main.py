import scipy
import numpy as np
import os.path

from ntpath import basename as basename

import astropy
from astropy.io import fits
from photutils.detection import DAOStarFinder,IRAFStarFinder

from ntpath import basename as basename

from matplotlib import pyplot as plt
from sys import argv


#Author-defined pckgs
import profileFitting as pffit

"""
OVERALL TODOS:
  - Some of the parameters in starfind should be made to update based on the measures of central tendancy above
  - Adapt source code to work for all sources within some brightness threshold
  - Build in support for changing plate scale (0.0317) for focal lengths
      config file?
  - Find a  way to extract bounds from the FITS file or make reasonable guesses
   Roy suggests supplying the function with bounds based on FITS data/physical constraints of our detectors
"""

inputPath=None
hdul=None

#function to load in fits data
def loadFits(inputPath):
  return fits.open(f'{inputPath}')[0]

#takes in file path or directory
if __name__ == '__main__': 
  try:
    inputPath = argv[1]
  except IndexError:
    print("No file path provided as argument to python script.\nExiting")
    exit(0)



  
  #empty list, to be filled with all inputed .fits files
  directoryFits=[]

  if inputPath.endswith(".fits") or inputPath.endswith(".fit"):
    #storing single file, finidng directory for that file
    directoryFits.append(loadFits(inputPath))
    inputPath=os.path.dirname(inputPath)
  else:
    #storing all fits to be processed
    for fitsFile in os.listdir(inputPath):
      if(fitsFile.endswith(".fits") or fitsFile.endswith(".fit")):
        directoryFits.append(fitsFile)


  #run calculations for every .fits file
  for fitsFile in directoryFits:
    ####Setting size of subFrames
    # leave as int!
    sfLength = 50

    #load in .fits file
    hdul=loadFits(inputPath+"/"+fitsFile)

    #Assumed to be in um
    pixSize = float(hdul.header['XPIXSZ'])

    data = hdul.data

    ####Find sources
    mean, median, std, max = np.mean(data), np.median(data), np.std(data), np.max(data)

    starFind = DAOStarFinder(threshold=median, fwhm=20.0, sky=500, exclude_border=True, brightest=10, peakmax=70000)
    sourcesList = starFind(data)

    #We'd need some logic here to determine which sources are good/bad to try to fit, for now just does brightest  

    for sourceID in [0]:

      ####Extract a subframe
      #Find center of first (brightest) source
      #  This code is not perfect, but works for a single source (here the 0th)
      xC, yC = sourcesList[0]['xcentroid'], sourcesList[0]['ycentroid']

      #Slicing (and matrices in general in python) are (row, col)
      subFrame = data[int(yC - sfLength):int(yC + sfLength),int(xC - sfLength):int(xC + sfLength)]

      ####Extract profile data (ideally centered on sources)     
      ##Extract Horizontal Data
      horizData = subFrame[sfLength,:]

      ##Extract Vertical Data
      vertiData = subFrame[:,sfLength]

      ##Extract Radial Profile Data
      # Assuming subFrame is square
      radialDataRaw = pffit.extract_radial_data(subFrame, xC=sfLength, yC=sfLength)[:sfLength]
      radialData = np.concatenate((radialDataRaw[::-1], radialDataRaw))

      ####Perform Fits
      xData = np.linspace(0, sfLength*2, sfLength*2)

      horizParams = pffit.fit_gaussian_1d(xData,horizData)
      vertiParams = pffit.fit_gaussian_1d(xData,vertiData)
      radialParams = pffit.fit_gaussian_1d(xData,radialData)
      
      #####Get Fitted Profiles
      horizFit  = pffit.gaussian_1d(xData, *horizParams)
      vertiFit  = pffit.gaussian_1d(xData, *vertiParams)
      radialFit = pffit.gaussian_1d(xData, *radialParams)

      ##Calculate Residuals
      ## SQRT( SUM( SQUARED_DIFFERENCES  ) ) / numDataPts
      horizResidual  = np.sqrt( sum( (horizFit  - horizData)  **2 ) ) / (2*sfLength)
      vertiResidual  = np.sqrt( sum( (vertiFit  - vertiData)  **2 ) ) / (2*sfLength)
      radialResidual = np.sqrt( sum( (radialFit - radialData) **2 ) ) / (2*sfLength)
      ####Convert the STD to FWHM and convert to " (arcsec) based on FITS header
      horizFWHM = 2.355*horizParams[1]  * 0.0317 * pixSize
      vertiFWHM = 2.355*vertiParams[1]  * 0.0317 * pixSize
      radFWHM   = 2.355*radialParams[1] * 0.0317 * pixSize
      
      print(f"Fits completed with the following residuals for {fitsFile}\nRadial: {radialResidual:0.3f}\nHorizontal: {horizResidual:0.3f}\nVertical: {vertiResidual:0.3f}\n")
      print(f"Radial FWHM: {radFWHM:0.3f}\nHorizontal FWHM: {horizFWHM:0.3f}\nVertical FWHM: {vertiFWHM:0.3f}")

      ####Generate Log
      with open(inputPath+"/FWHMscript-output-log.txt", "a") as f:
        #Add data block line break
        f.write("========================================\n")

        #Source ID (internal to this code), ordered by brightness of subframes
        f.write(f"Source ID: {sourceID}\n") # Come back to this when loopen

        #File Name
        f.write(f"File Name: {basename(fitsFile)} \n ")

        #Observation Time
        datePart, timePart = hdul.header['DATE-OBS'].split("T")
        year, month, day = datePart.split("-")
        hour, minute = timePart.split(":")[:2]
        obsTime=f"{year}{month}{day}T{hour}{minute}"
        
        f.write(f"Obs Start Time: {obsTime}\n")

        #Camera used (with pixel size)
        f.write(f"Camera Used (pixel size): {hdul.header['INSTRUME']} ({hdul.header['XPIXSZ']}(um))\n")

        #Per type of profile
        ##Write all HORIZONTAL fit parameters with residuals and FWHM
        f.write(f"Horizontal Fit\nmu:\t\t{horizParams[0]}\nsigma:\t\t{horizParams[1]}\namplitude:\t{horizParams[2]}\n"+\
                f"offset:\t\t{horizParams[3]}\nResidual:\t{horizResidual}\nFWHM:\t\t{horizFWHM}\n")
        ##Write all VERTICAL fit parameters with residuals and FWHM
        f.write(f"Vertical Fit\nmu:\t\t{vertiParams[0]}\nsigma:\t\t{vertiParams[1]}\namplitude:\t{vertiParams[2]}\n"+\
                f"offset:\t\t{vertiParams[3]}\nResidual:\t{vertiResidual}\nFWHM:\t\t{vertiFWHM}\n")
        ##Write all RADIAL fit parameters with residuals and FWHM
        f.write(f"Radial Fit\nmu:\t\t{radialParams[0]}\nsigma:\t\t{radialParams[1]}\namplitude:\t{radialParams[2]}\n"+\
                f"offset:\t\t{radialParams[3]}\nResidual:\t{radialResidual}\nFWHM:\t\t{radFWHM}\n")


      ####Generate Plots
      fig,charts =plt.subplots(2,2, figsize=(10,8))

      #display horizontal
      charts[0,0].plot(xData, horizData, 'ko', markersize=2)
      charts[0,0].plot(xData, horizFit,'tab:blue', linestyle='dashed',)
      charts[0,0].set_title("Horizontal Fit")
      charts[0,0].set(xlabel='Pixel in SubFrame',ylabel='Counts')
      charts[0,0].grid(1)

      #display vertical
      charts[0,1].plot(xData, vertiData, 'ko', markersize=2)
      charts[0,1].plot(xData, vertiFit,'tab:red', linestyle='dashed',)
      charts[0,1].set_title("Vertical Fit")
      charts[0,1].set(xlabel='Pixel in SubFrame',ylabel='Counts')
      charts[0,1].grid(1)

      #display radial
      charts[1,0].plot(xData,radialData, 'ko', markersize=2)
      charts[1,0].plot(xData,radialFit,'tab:purple', linestyle='dashed')
      charts[1,0].set_title("Radial Fit")
      charts[1,0].set(xlabel='Pixel in SubFrame',ylabel='Counts')
      charts[1,0].grid(1)

      #display subframe
      charts[1,1].imshow(subFrame, cmap='gray')
      charts[1,1].plot(sfLength,sfLength, 'rx')
      charts[1,1].set_title("Source SubFrame")

      plt.suptitle(f"FWHM Curve Fitting for Source ID: {sourceID}\n{fitsFile}")
      plt.tight_layout()
      plt.savefig("{}/{}_{}.png".format(inputPath,fitsFile[:-5], sourceID))