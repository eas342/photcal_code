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
    outList = vquery.query_region(coor,radius=5*u.arcsec,catalog="I/345/gaia2")
    outDat = outList[0]
    
    for oneColumn in dat.colnames:
        outDat["{}_PS".format(oneColumn)] = dat[oneColumn]
    
    for oneFormat in ['fits','csv']:
        outDat.write('lists/gaia_coord/gaia_{}.{}'.format(basePrefix,oneFormat),
                overwrite=True)
    
    