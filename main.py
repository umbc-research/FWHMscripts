import scipy
import numpy as np
import os.path
import datetime as datetime
import csv
from ntpath import basename as basename

import astropy
from astropy.io import fits
from photutils.detection import DAOStarFinder,IRAFStarFinder

from ntpath import basename as basename

from matplotlib import pyplot as plt
from sys import argv


#Author-defined pckgs
import profileFitting as pffit

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




  #generate log file for current run
  fields=['file name','source id','obs time', 'camera used (pixel size)',\
        'horiztonal mu', 'horiztonal sigma', 'horiztonal amplitude','horiztonal offset',\
        'horizontal R2', 'horizontal FWHM pixel','horizontal FWHM arcsecond',\
        'vertical mu', 'vertical  sigma', 'vertical  amplitude','vertical  offset',\
        'vertical  R2', 'vertical  FWHM pixel','vertical  FWHM arcsecond',\
        'radial mu', 'radial sigma', 'radial amplitude','radial offset',\
        'radial R2', 'radial FWHM pixel','radial FWHM arcsecond'\
                ]


  dTime=f'{datetime.datetime.now()}'
  year=f'{dTime[:4]}'
  month=f'{dTime[5:7]}'
  day=f'{dTime[8:10]}'
  hour=f'{dTime[11:13]}'
  minute=f'{dTime[14:16]}'
  second=f'{dTime[17:19]}'
  runTime =f'{year}{month}{day}T{hour}{minute}{second}'
  with open(f'{inputPath}/FWHMscript-output-log-{runTime}.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(fields)


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

      ##Calculate R^2
      horizResList = (horizFit  - horizData)
      horizSumSqrsRes=np.sum(horizResList*horizResList)
      horizTotSumSqrs=np.sum((horizData-np.mean(horizData))**2)
      horizR2=1.0-(horizSumSqrsRes/horizTotSumSqrs)

      vertiResList = (vertiFit  - vertiData)
      vertiSumSqrsRes=np.sum(vertiResList*vertiResList)
      vertiTotSumSqrs=np.sum((vertiData-np.mean(vertiData))**2)
      vertiR2=1.0-(vertiSumSqrsRes/vertiTotSumSqrs)

      radialResList = (radialFit  - radialData)
      radialSumSqrsRes=np.sum(radialResList*radialResList)
      radialTotSumSqrs=np.sum((radialData-np.mean(radialData))**2)
      radialR2=1.0-(radialSumSqrsRes/radialTotSumSqrs)

      # residuals code, no longer used
      ## SQRT( SUM( SQUARED_DIFFERENCES  ) ) / numDataPts
      # horizResidual  = np.sqrt( sum( (horizFit  - horizData)  **2 ) ) / (2*sfLength)
      # vertiResidual  = np.sqrt( sum( (vertiFit  - vertiData)  **2 ) ) / (2*sfLength)
      # radialResidual = np.sqrt( sum( (radialFit - radialData) **2 ) ) / (2*sfLength)

      
      ####Convert the STD to FWHM and convert to " (arcsec) based on FITS header
      horizFWHMpix=2.355*horizParams[1]
      vertiFWHMpix=2.355*vertiParams[1]
      radialFWHMpix=2.355*radialParams[1]

      horizFWHMarc =  horizFWHMpix * 0.0317 * pixSize
      vertiFWHMarc =  vertiFWHMpix  * 0.0317 * pixSize
      radialFWHMarc = radialFWHMpix* 0.0317 * pixSize
      
      print(f"Fits completed with the following R^2 s for {fitsFile}\nHorizontal: {horizR2:0.3f}\nVertical: {vertiR2:0.3f}\nRadial: {radialR2:0.3f}\n")
      print(f"Horizontal FWHM(arcseconds): {horizFWHMarc:0.3f}\nVertical FWHM(arcseconds): {vertiFWHMarc:0.3f}\nRadial FWHM (arcseconds): {radialFWHMarc:0.3f}\n")

      ####Generate Log
      uselesspart, data = f"{hdul.header['SIMPLE*']}".split("FITS: ")
      timeInfo, uselesspart2 = data.split("E")
       
      datePart=timeInfo[:10]
      timePart=timeInfo[11:].strip()
      month, day, year = datePart.split("/")
      hour, minute,second = timePart.split(":")
      obsTime=f"{year}{month}{day}T{hour}{minute}{second}"
        
      data=[f'{basename(fitsFile)}',f'{sourceID}',f'{obsTime}',f'{hdul.header["INSTRUME"]} ({hdul.header["XPIXSZ"]} (um))',\
        f'{horizParams[0]}', f'{horizParams[1]}',f'{horizParams[2]}',f'{horizParams[3]}',\
        f'{horizR2}', f'{horizFWHMarc}', f'{horizFWHMpix}', \
        f'{vertiParams[0]}', f'{vertiParams[1]}',f'{vertiParams[2]}',f'{vertiParams[3]}',\
        f'{vertiR2}', f'{vertiFWHMarc}', f'{vertiFWHMpix}', \
        f'{radialParams[0]}', f'{radialParams[1]}',f'{radialParams[2]}',f'{radialParams[3]}',\
        f'{radialR2}', f'{radialFWHMarc}', f'{radialFWHMpix}']
      
      with open(f'{inputPath}/FWHMscript-output-log-{runTime}.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(data)
      
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
      plt.savefig("{}/{}_{}_{}.png".format(inputPath,fitsFile[:-5], sourceID,runTime))
  