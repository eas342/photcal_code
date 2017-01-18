from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
import es_gen
from astropy.convolution import Gaussian1DKernel, convolve
import numpy as np

class clusterSpec(object):
    def __init__(self,kernelWidth=2.5,
                redFile='../hectospec_data/2016.0207/reduction/0102/spHect-ngc2420_cat_1.2468-0102.fits'):
        hdulist = fits.open(redFile)
        self.fibinfo = Table(hdulist[5].data)
        self.makeCleanNames()
        self.goodfib = np.array(es_gen.es_strmatch('O-*',self.fibinfo['OBJTYPE']))
        self.wave2D = hdulist[0].data
        self.flux2D = hdulist[1].data
        self.fluxErr2D = 1/np.sqrt(hdulist[2].data)
        self.kernelWidth=kernelWidth
    
    def getbyObjName(self,objName):
        """ Search by object name
        """
        foundPts = self.fibinfo['cleanName'] == objName
        nFound = np.sum(foundPts)
        if np.sum(foundPts) == 1:
            return self.fibinfo[foundPts]
        else:
            print("Found "+str(nFound)+" objects for "+objName)
    
    def makeCleanNames(self):
        """
        Makes clean versions of the objtypes (no plus signs for use in file names)
        """
        cleanNames = []
        problemStrings = ['O-','-','+']
        replaceStrings = ['O_','_m','_p']
        for oneRow in self.fibinfo:
            newString = oneRow['OBJTYPE']
            for problemString, replaceString in zip(problemStrings,replaceStrings):
                newString = replaceString.join(newString.split(problemString))
            cleanNames.append(newString)
        self.fibinfo['cleanName'] = cleanNames
        
    def prepForMK(self):
        """
        Prepares the spectra for MKClass
        """
        datkern = Gaussian1DKernel(stddev=self.kernelWidth)
        for currentFib in self.goodfib:
            #fig, ax = plt.subplots(figsize=(15,4))
            x = self.wave2D[currentFib,:]
            y = self.flux2D[currentFib,:]
            yconv = convolve(y, datkern, boundary='extend')
            ynorm = yconv / np.median(yconv)
            keepPts = (x > 3800) & (x < 5600)
            outX = x[keepPts]
            #CorrectionFunc = np.interp(outX,xinterp,ymod)
            outY = ynorm[keepPts]
            #outY = ynorm[keepPts] / CorrectionFunc
            #plt.plot(outX,outY)
    
            outTable = Table()
            outTable['Wavelength_A']=outX
            outTable['Flux']=outY
            objFileName = self.fibinfo['cleanName'][currentFib]
            ascii.write(outTable,'output_flux_cal/'+objFileName+'.txt',format='fixed_width_no_header',delimiter=' ')
        