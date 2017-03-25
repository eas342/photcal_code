from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, vstack
import es_gen
from astropy.convolution import Gaussian1DKernel, convolve
import numpy as np
import yaml
import multi_module
import matplotlib.pyplot as plt
import pdb
import os

## Remove the warnings on the robust polynomial fitting
import warnings
warnings.simplefilter('ignore', np.RankWarning)

class clusterSpec(object):
    def __init__(self,kernelWidth=2.5,src='NGC 2420',indices=[0]):
        self.src = src
        self.srcFileName = src.replace(" ","_")
        self.redFiles = []
        for oneInd in indices:
            redFile = multi_module.getRedFile(src,dataType='hectoData',index=oneInd)
            self.redFiles.append(redFile)
            
            thisFibinfo, thisWave2D, thisFlux2D, thisFluxErr2D = self.getFITSData(redFile)
            if oneInd >= 1:
                fibinfo = vstack([fibinfo,thisFibinfo])
                wave2D = np.vstack([wave2D,thisWave2D])
                flux2D = np.vstack([flux2D,thisFlux2D])
                fluxErr2D = np.vstack([fluxErr2D,thisFluxErr2D])
            else:
                fibinfo, wave2D, flux2D, fluxErr2D = thisFibinfo, thisWave2D, thisFlux2D, thisFluxErr2D
            
        self.fibinfo = fibinfo
        self.wave2D = wave2D
        self.flux2D = flux2D
        self.fluxErr2D = fluxErr2D
        
        self.goodfib = np.array(es_gen.es_strmatch('O-*',self.fibinfo['OBJTYPE']))
        self.makeCleanNames()
        self.kernelWidth=kernelWidth
    
    def getFITSData(self,fileName):
        """ 
            Gets Hectospec fiber info, wavelength flux and error
        """
        hdulist = fits.open(fileName)
        fibinfo = Table(hdulist[5].data)
        wave2D = hdulist[0].data
        flux2D = hdulist[1].data
        fluxErr2D = 1/np.sqrt(hdulist[2].data)
        return fibinfo, wave2D, flux2D, fluxErr2D
    
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
        problemStrings = ['O-','-','+',r' ']
        replaceStrings = ['O_','_m','_p','']
        for oneRow in self.fibinfo:
            newString = oneRow['OBJTYPE']
            for problemString, replaceString in zip(problemStrings,replaceStrings):
                newString = replaceString.join(newString.split(problemString))
            
            cleanNames.append(newString)
        self.fibinfo['cleanName'] = cleanNames
        
    def prepForMK(self,showSteps=False,nset = 3,waveStart=3800,waveEnd=5572):
        """
        Prepares the spectra for MKClass
        
        Parameters
        -----------------
        showSteps: bool
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
            targFileName = self.fibinfo['cleanName'][currentFib]
            
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
            
            if showSteps == True:
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
                
                fig.savefig('plots/spec_rect_process/'+targFileName+'_spec.pdf')
                fig.show()
                pdb.set_trace()
            
            
            ## Choose the directories in which to save plots
            plotSaveDir = os.path.join('plots','spec_rect',self.srcFileName)
            rectSpecSaveDir = os.path.join('output_rectified',self.srcFileName)
            saveDirs = [plotSaveDir,rectSpecSaveDir]
                        
            for oneDir in saveDirs:
                if os.path.exists(oneDir) == False:
                    os.mkdir(oneDir)
            
            
            fig,ax = plt.subplots(figsize=(8,4))
            ax.plot(outX,ynorm2)
            ax.set_xlabel('Wavelength ($\AA$)')
            ax.set_ylabel('Normalized Flux')
            fig.savefig(os.path.join(plotSaveDir,targFileName+'_spec.pdf'))
            
            t = Table()
            t['Wavelength'] = outX
            t['Normalized Flux'] = ynorm2
            ascii.write(t,os.path.join(rectSpecSaveDir,targFileName+'_spec.txt'),format='fixed_width_no_header',
                        delimiter=' ',overwrite=True)
            
            plt.close('all')
            
        