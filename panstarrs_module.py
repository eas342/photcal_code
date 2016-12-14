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

if __name__ == "__main__":
    cP = clusterPhot()
    cS = hm.clusterSpec()
    