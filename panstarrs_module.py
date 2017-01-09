from astropy.io import ascii
from astropy.io import fits
import matplotlib.pyplot as plt
import hecto_module as hm
import numpy as np

class clusterPhot(object):
    """ Creates a cluster object for Pan-Starrs Photometry """
    
    def __init__(self,photFile='../pan_starrs/NGC2420.txt'):
        self.photFile = photFile
        self.dat = ascii.read(photFile)
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
    
    def get_cluster_pt(self,dist=5,ra='ra',dec='dec',racen=114.5958,deccen=21.573):
        """ 
        Gets the cluster points from the distance to center
        
        Parameters
        --------------
        dist: float
            Distance in arc-minutes
        """
        deltaRA = racen - self.dat[ra]
        deltaDEC = deccen - self.dat[dec]
        deltaDistApprox = np.sqrt(deltaRA**2 / np.sin(self.dat[dec])**2 + deltaDEC**2)
        self.cpoints = deltaDistApprox < dist/60.
        
    
    def plot_cm(self,color1='g',color2='r',mag='g',):
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
        fig, ax = plt.subplots()
        cdat = self.dat[self.cpoints]
        
        ax.plot(cdat[color1] - cdat[color2],cdat[mag],'.',rasterized=True)
        ax.set_xlabel(color1+' - '+color2)
        ax.set_ylabel(mag)
        ax.invert_yaxis()
        
        self.ax = ax
        self.fig = fig
        #fig.show()

if __name__ == "__main__":
    cP = clusterPhot()
    cS = hm.clusterSpec()
    