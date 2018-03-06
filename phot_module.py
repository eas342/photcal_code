from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import hecto_module as hm
import numpy as np
import multi_module
import pdb
import os


## Dictionary to get correct coordinate column coordinates
coordinateDict = {'_RA':'ra','Ra (deg)':'ra',
                  '_DE':'dec','Dec (deg)':'dec'}



class clusterPhot(object):
    """ Creates a cluster object for Pan-Starrs Photometry
    Also works with photometry for other sources
    Specific the photometry with PhotType
    
    Paraemters
    ---------------
    photType: str
        Photometry type
         - `panStarrsData` - Pan Starrs griz
         - `UKIRTData` - UKIRT JHK
    """
    
    def __init__(self,src='NGC 2420',photType='panStarrsData'):
        self.src = src
        photFile = multi_module.getRedFile(src,dataType=photType)
        racen = multi_module.getClusterInfo(src,'RA') ## Decimal degrees
        deccen = multi_module.getClusterInfo(src,'Dec') ## Decimal degrees
        self.photFile = photFile
        if os.path.splitext(self.photFile)[-1] == '.fits':
            HDUList = fits.open(self.photFile)
            self.dat = Table(HDUList[1].data)
            for oneKey in coordinateDict.keys():
                if oneKey in self.dat.colnames:
                    self.dat[coordinateDict[oneKey]] = self.dat[oneKey]
            
            HDUList.close()
        else:
            self.dat = ascii.read(self.photFile)
        
        self.photCoor = SkyCoord(ra=self.dat['ra'] * u.degree,
                                 dec=self.dat['dec'] * u.degree)
        
        self.racen = racen ## Decimal degrees
        self.deccen = deccen ## Decimal degrees
        self.cluster_cen = SkyCoord(ra=racen * u.degree,dec=deccen * u.degree)
        
        self.clusterRad = multi_module.getClusterInfo(src,'Cluster Rad (arcmin)') * u.arcmin
        self.get_cluster_pt()
        
        ## Threshold for target matching
        self.targMatch = 0.2 * u.arcsec
    
    def lookup_src(self,ra,dec):
        """ Look up a source from the RA and Dec. Expects these to both be in degrees"""
        
        coor = SkyCoord(ra=ra * u.degree,dec=dec * u.degree)
        rowIndex, minimumDist, d2d = coor.match_to_catalog_sky(self.photCoor)
        
        if minimumDist > self.targMatch:
            print('Warning, closest source is '+str(minimumDist.arcsec)+' arcsec for'+str(ra)+' '+str(dec))
            return None
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
        dist = self.cluster_cen.separation(self.photCoor)
        
        self.cpoints = dist < self.clusterRad
        
    
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
    