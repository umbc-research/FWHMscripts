# FWHM Scripts

Authors:Olivia Chiarini and Roy Prouty
  
  
Research Advisors: Dr. Meyer and Roy Prouty

  
Purpose: To analyze FWHM data from .fits files given by UMBC's Obsertvatory, compare these to SHARPCAP to see how they are getting their histograms for FWHM and focusing analyizse. Future goals are to generalize this to be used for data processing at the observatory as a whole.
  

OVERALL TODOS:

- Some of the parameters in starfind should be made to update based on the measures of central tendancy above

- Adapt source code to work for all sources within some brightness threshold

- Build in support for changing plate scale (0.0317) for focal lengths

    config file?

- Find a way to extract bounds from the FITS file or make reasonable guesses

    Roy suggests supplying the function with bounds based on FITS data/physical constraints of our detectors