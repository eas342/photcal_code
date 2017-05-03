from astropy.io import ascii
from astropy.io import fits
import matplotlib.pyplot as plt
import hecto_module as hm
import numpy as np
import multi_module
import pdb

class clusterPhot(object):
    """ Creates a cluster object for Pan-Starrs Photometry """
    
    def __init__(self,src='NGC 2420',clusterRad=7.5):
        self.src = src
        photFile = multi_module.getRedFile(src,dataType='panStarrsData')
        racen = multi_module.getClusterInfo(src,'RA') ## Decimal degrees
        deccen = multi_module.getClusterInfo(src,'Dec') ## Decimal degrees
        self.photFile = photFile
        self.dat = ascii.read(photFile)
        self.racen = racen ## Decimal degrees
        self.deccen = deccen ## Decimal degrees
        self.clusterRad = multi_module.getClusterInfo(src,'Cluster Rad (arcmin)')
        self.get_cluster_pt()
    
    def lookup_src(self,ra,dec):
        """ Look up a source from the RA and Dec. Expects these to both be in degrees"""
        deltaRA = ra - self.dat['ra']
        deltaDEC = dec - self.dat['dec']
        deltaDistApprox = np.sqrt(deltaRA**2 / np.sin(self.dat['dec'])**2 + deltaDEC**2)
        rowIndex, minimumDist = np.argmin(deltaDistApprox), np.min(deltaDistApprox)
        if minimumDist > 0.08 / 3600.:
            print('Warning, closest source is '+str(minimumDist*3600.)+' arcsec for'+str(ra)+' '+str(dec))
        else:
            return self.dat[rowIndex]
    
    def get_cluster_pt(self,ra='ra',dec='dec'):
        """ 
        Gets the cluster points from the distance to center
        
        Parameters
        --------------
        dist: float
            Distance in arc-minutes
        """
        deltaRA = self.racen - self.dat[ra]
        deltaDEC = self.deccen - self.dat[dec]
        deltaDistApprox = np.sqrt(deltaRA**2 * np.cos(self.dat[dec])**2 + deltaDEC**2)
        dist= self.clusterRad / 60.
        self.cpoints = deltaDistApprox < dist
        
    
    def plot_cm(self,color1='g',color2='r',mag='g',figsize=None):
        """ 
        Plots a color-magnitude diagram from the photometry 
        
        Parameters
        --------------
        color1: str
            First photometric band
        color2: str
            Second photometric band
        mag: str
            Magnitude photometric band
        """
        
        fig, ax = plt.subplots(figsize=figsize)
        cdat = self.dat[self.cpoints]
        
        ax.plot(cdat[color1] - cdat[color2],cdat[mag],'.',rasterized=True,label='')
        ax.set_xlabel(color1+' - '+color2)
        ax.set_ylabel(mag)
        ax.invert_yaxis()
        
        self.ax = ax
        self.fig = fig
    
    def plot_fov(self):
        """
        Plots the cluster stars in RA and DEC
        """
        fig, ax = plt.subplots()
        cdat = self.dat[self.cpoints]
        
        ax.plot(cdat['ra'],cdat['dec'],'.',rasterized=True,label='')
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        
        ### Show the NIRCam FOV
        NIRCamFOV = 2.16 / 60. ## deg
        deltaX = NIRCamFOV * 0.5/ np.cos(self.deccen * np.pi/180.) ## Distances in RA/Dec, deg
        deltaY = NIRCamFOV * 0.5 ## Distances in RA/dec deg

        coordA = self.racen - 0.02, self.deccen + 0.03
        coordB = coordA[0], coordA[1] - NIRCamFOV - 42/3600.
        
        for ncoord, FOVname in zip([coordA,coordB],['A','B']):
            xpt = ncoord[0] + np.array([-deltaX,deltaX,deltaX,-deltaX,-deltaX])
            ypt = ncoord[1] + np.array([-deltaY,-deltaY,deltaY,deltaY,-deltaY])
            ax.plot(xpt,ypt,linewidth=1,label='Module '+FOVname)
        #ax.legend(frameon=False,loc='best')
        
        self.ax = ax
        self.fig = fig

if __name__ == "__main__":
    cP = clusterPhot()
    cS = hm.clusterSpec()
    