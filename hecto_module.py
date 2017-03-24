from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
import es_gen
from astropy.convolution import Gaussian1DKernel, convolve
import numpy as np
import yaml
import multi_module
import matplotlib.pyplot as plt
import pdb

## Remove the warnings on the robust polynomial fitting
import warnings
warnings.simplefilter('ignore', np.RankWarning)

class clusterSpec(object):
    def __init__(self,kernelWidth=2.5,src='NGC 2420'):
        self.src = src
        self.redFile = multi_module.getRedFile(src,dataType='hectoData')
        hdulist = fits.open(self.redFile)
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
        
    def prepForMK(self,showFig=False,nset = 3,waveStart=3800,waveEnd=5600):
        """
        Prepares the spectra for MKClass
        
        Parameters
        -----------------
        showFig: bool
            Show a figure for the spectra rectification
        nset: int
            Number of regions over which to do polynomial fitting
        
        """
        datkern = Gaussian1DKernel(stddev=self.kernelWidth)
        
        ## Wavelength locations for the regions over which to do polynomials
        waveLocs = np.linspace(waveStart,waveEnd,nset+1)
        
        for currentFib in self.goodfib:
            #fig, ax = plt.subplots(figsize=(15,4))
            x = self.wave2D[currentFib,:]
            y = self.flux2D[currentFib,:]
            yconv = convolve(y, datkern, boundary='extend')
            
            ## Points with useful features for classification by mkClass
            keepPts = (x > waveStart) & (x < waveEnd)
            outX = x[keepPts]
            trimY = yconv[keepPts]
            
            ## Do a fit in log space so it won't give you zeros
            yLog = np.log10(trimY)
            
            modelF1 = np.zeros_like(outX)
            
            setPtsArr = []
            for indSet in range(nset):
                setPts = (outX > waveLocs[indSet]) & (outX <= waveLocs[indSet+1])
                setPtsArr.append(setPts)
                setPoly = es_gen.robust_poly(outX[setPts],yLog[setPts],8,sigreject=2.)
                modelF1[setPts] = 10**np.polyval(setPoly,outX[setPts])
            
            ynorm1 = trimY / modelF1
            
            ## Choose the points that are above 1.0 for the continuum
            contPoints = (ynorm1 > 1.0)
            
            modelF2 = np.zeros_like(outX)
            fitPtsArr = []
            for indSet in range(nset):
                setPts = (setPtsArr[indSet])
                fitPts = setPts & contPoints ## continuum and part of the set
                fitPtsArr.append(fitPts)
                ## add in the point previous point for continuity
                if indSet >= 1:
                    ## Get the maximum index of the previous set of fit points
                    prevPoint = np.argmax(outX[fitPtsArr[indSet-1]])
                    
                setPoly = es_gen.robust_poly(outX[fitPts],ynorm1[fitPts],8,sigreject=3.)
                modelF2[setPts] = np.polyval(setPoly,outX[setPts])
            ynorm2 = ynorm1 / modelF2
            
            if showFig == True:
                plt.close('all')
                fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
                ax1.semilogy(outX,trimY)
                ax1.plot(outX,modelF1)
                ax1.set_ylabel('Flux')
                ax2.plot(outX,ynorm1)
                ax2.plot(outX,modelF2)
                ax2.set_ylabel('Norm Flux')
                ax3.plot(outX,ynorm2)
                ax3.set_ylabel('Norm Flux')
                ax3.set_xlabel('Wavelength ($\AA$)')
                
                fig.savefig('plots/spec_rect/'+self.fibinfo['cleanName'][currentFib]+'_spec.pdf')
                fig.show()
            pdb.set_trace()

            #CorrectionFunc = np.interp(outX,xinterp,ymod)
            outY = ynorm[keepPts]
            #outY = ynorm[keepPts] / CorrectionFunc
            #plt.plot(outX,outY)
    
            outTable = Table()
            outTable['Wavelength_A']=outX
            outTable['Flux']=outY
            objFileName = self.fibinfo['cleanName'][currentFib]
            ascii.write(outTable,'output_flux_cal/'+objFileName+'.txt',format='fixed_width_no_header',delimiter=' ')
        