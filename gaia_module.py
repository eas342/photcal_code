"""
## This module is designed to retrieve Gaia positions, parallaxes
proper motions and magnitudes
"""

from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier 
import astropy.units as u 
from astropy.time import Time
from astropy.io import ascii, fits
from distutils.version import LooseVersion
import logging
import pdb
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

#if LooseVersion(astropy.__version__) < LooseVersion("3.0"):
#    logging.error("Need Astropy >=3.0 for this script. You'll also need Python >3")


def get_gaia():
    """ Gets the Gaia info for a table """
    
    basePrefix = "lris_targsNGC2506"
    fullDat = ascii.read('../pan_starrs/pro/output/{}.csv'.format(basePrefix))
    pts = fullDat['GROUP'] == 1
    dat = fullDat[pts]
    
    coor = SkyCoord(dat['RA'],dat['DEC'],unit=(u.deg,u.deg))
    
    ## Set up Vizier Query
    vquery = Vizier(row_limit = 1)
    outList = vquery.query_region(coor,radius=1*u.arcsec,catalog="I/345/gaia2")
    outDat = outList[0]
    
    for oneColumn in dat.colnames:
        outDat["{}_PS".format(oneColumn)] = dat[oneColumn]
    
    for oneFormat in ['fits','csv']:
        outDat.write('lists/gaia_coord/gaia_{}.{}'.format(basePrefix,oneFormat),
                overwrite=True)
    
# def show_offsets():
#     """ Show the offsets between Pan-Starrs and Gaia DR2 """
#     dat = Table.read('lists/gaia_coord/gaia_lris_targsNGC2506.fits')
#     diffRA = (dat['RA_ICRS'] - dat['RA_PS']) * 3600.
#     diffDec = (dat['DE_ICRS'] - dat['DEC_PS']) * 3600.
#     plt.plot(diffRA,diffDec,'.')
#     plt.xlabel('')
#     plt.ylabel()
#     plt.show()

def compare_ps(compareParameter='mag'):
    """ 
    Compare Gaia DR2 and PS
    
    Parameters
    ----------
    compareParameter: str
        If 'mag', show the g - G magnitude difference
        If 'dist', show the offsets in coordinates between Pan-Starrs
        and Gaia DR2
    """
    
    dat = Table.read('lists/gaia_coord/gaia_lris_targsNGC2506.fits')
    
    diffRA = (dat['RA_ICRS'] - dat['RA_PS']) * 3600.
    diffDec = (dat['DE_ICRS'] - dat['DEC_PS']) * 3600.
    
    closePt = np.sqrt(diffRA**2 + diffDec**2) < 0.12
    
    color = dat['G_PS'] - dat['R_PS']
    deltaMag = dat['G_PS'] - dat['Gmag']
    
    if compareParameter == 'mag':
        x, y = color, deltaMag
        xLabel = 'g$_{PS}$ - r$_{PS}$'
        yLabel = 'g$_{PS}$ - G$_{Gaia}$'
    else:
        x, y = diffRA, diffDec
        xLabel = '$\Delta$ RA (arcsec)'
        yLabel = '$\Delta$ Dec (arcsec)'
    
    for ind, oneType in enumerate(['Close','Far']):
        if oneType == 'Close':
            pts = closePt
        else:
            pts = (closePt == False)
        plt.plot(x[pts],y[pts],'.',label=oneType)
    
    plt.legend()
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
