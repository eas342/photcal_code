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
    
    outDat['diffRA'] = (outDat['RA_ICRS'] - outDat['RA_PS']) * 3600.
    outDat['diffDEC'] = (outDat['DE_ICRS'] - outDat['DEC_PS']) * 3600.
    
    for oneFormat in ['fits','csv']:
        outDat.write('lists/gaia_coord/gaia_{}.{}'.format(basePrefix,oneFormat),
                overwrite=True)
    

def group_sources(dat,groupParameter='ra'):
    """ 
    Group the cluster size
    
    Parameters
    -----------
    groupParameter: str
        If 'all', group by several criteria
    """
    if groupParameter == 'dist':
        closePt = np.sqrt(dat['diffRA']**2 + dat['diffDEC']**2) < 0.12
    elif groupParameter == 'ra':
        closePt = dat['diffRA']  < 0.05
    elif groupParameter == 'all':
        closePt = np.ones(len(dat),dtype=np.bool)
        for oneParam in ['pmRA','pmDE']:
            med = np.median(dat[oneParam])
            mad = np.median(np.abs(dat[oneParam] - med))
            diff = np.abs(dat[oneParam] - med)
            closePt = closePt & (diff < 5.0 * mad)
        
    else:
        print('Unrecognized grouping parameter')
        closePt = np.isfinite(diffRA)
    
    return closePt
    

def compare_ps(compareParameter='mag',
               groupParameter='all'):
    """ 
    Compare Gaia DR2 and PS
    
    Parameters
    ----------
    compareParameter: str
        If 'mag', show the g - G magnitude difference
        If 'dist', show the offsets in coordinates between Pan-Starrs
        and Gaia DR2
        If 'pm', show the proper motions
    groupParameter: str
        If 'dist', group the sources by Gaia-PS distance
        If 'ra', group by RA
    """
    
    dat = Table.read('lists/gaia_coord/gaia_lris_targsNGC2506.fits')

    
    closePt = group_sources(dat,groupParameter=groupParameter)
    
    color = dat['G_PS'] - dat['R_PS']
    deltaMag = dat['G_PS'] - dat['Gmag']
    
    if compareParameter == 'mag':
        x, y = color, deltaMag
        xLabel = 'g$_{PS}$ - r$_{PS}$'
        yLabel = 'g$_{PS}$ - G$_{Gaia}$'
    elif compareParameter == 'dist':
        x, y = dat['diffRA'], dat['diffDEC']
        xLabel = '$\Delta$ RA (arcsec)'
        yLabel = '$\Delta$ Dec (arcsec)'
    elif compareParameter == 'pm':
        x, y = dat['pmRA'], dat['pmDE']
        xLabel = 'Proper Motion RA (mas/yr)'
        yLabel = 'Proper Motion Dec (mas/yr)'
    
    for ind, oneType in enumerate(['Cluster','Excluded']):
        if oneType == 'Cluster':
            pts = closePt
        else:
            pts = (closePt == False)
        plt.plot(x[pts],y[pts],'.',label=oneType)
    
    plt.legend()
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
def make_autoslit():
    """ Make a file for autoslit """
    