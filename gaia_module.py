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
import astropy.table
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

#import np.core.defchararray.add as stringadd

#if LooseVersion(astropy.__version__) < LooseVersion("3.0"):
#    logging.error("Need Astropy >=3.0 for this script. You'll also need Python >3")


def get_gaia(basePrefix="lris_targsNGC2506",
             groupInclude=[1]):
    """ Gets the Gaia info for a table 
    Parameters
    --------------
    basePrefix: str
        "lris_targsNGC2506" is the list of solar analog candidates
        "lris_alignmentNGC2506" is the list of alignment stars
    groupInclude: list
        List of integers describing the groups of stars to include
        These are priority levels on stars with 1 the highest
    """
    
    
    fullDat = ascii.read('../pan_starrs/pro/output/{}.csv'.format(basePrefix))
    pts = np.zeros(len(fullDat),dtype=np.bool)
    for oneGroup in groupInclude:
        pts = pts | (fullDat['GROUP'] == oneGroup)
    
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
    

def group_sources(dat,groupParameter='all'):
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
    

def cluster_distance(basePrefix="full_targ_lists_NGC2506"):
    """
    Get the Gaia distance to the cluster from the main sequence stars
    """
    msDat = Table.read("lists/gaia_coord/gaia_{}.{}".format(basePrefix,"fits"))
    goodp = msDat['e_Plx'] < 0.05
    
    medianPlx = np.nanmedian(msDat['Plx'][goodp])
    errPlx = scipy.stats.sem(msDat['Plx'][goodp])
    
    print("Plx= {} +/- {}".format(medianPlx,errPlx))
    dist = 1000./medianPlx
    errDist = dist * errPlx/medianPlx
    print("Dist= {} pc +/- {}".format(dist,errDist))
    
    distMod0 = 5. * np.log10(dist/10.)
    errDistMod0 = 5. * errDist / (np.log(10.) * dist)
    
    print("Dist Mod = {} +/- {}".format(distMod0,errDistMod0))
    

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
    
def make_autoslit(runRound=1):
    """ Make a file for autoslit
    
    Parameters
    ----------
    runRound: int
        Round to be run
     """
    
    ## read in the file
    fullDat = Table.read('lists/gaia_coord/gaia_lris_targsNGC2506.fits')
    
    ## only include cluster points
    clusterPt = group_sources(fullDat,groupParameter='all')
    dat = fullDat[clusterPt]
    
    ## Read in 6 alignment stars
    #alignDat = Table.read('lists/gaia_coord/subset_1_alignment_stars.csv')
    alignDatFull = Table.read('lists/gaia_coord/gaia_lris_alignmentNGC2506.fits')
    goodAlign = np.isfinite(alignDatFull['pmRA']) & np.isfinite(alignDatFull['pmDE'])
    alignDat = alignDatFull[goodAlign]
    
    ## get the alignment
    alignName = []
    for ind,oneRow in enumerate(alignDat):
        alignName.append("Align-{:02d}-g={:.2f}".format(ind+1,oneRow['G_PS']))
    
    ## Read this from FILE if generalizing to other clusters
    grSolarColor = 0.43
    
    t = Table()
    tAlign = Table()
    
    ## Formerly, give them names from coordinates
    coor = SkyCoord(dat['RA_ICRS'],dat['DE_ICRS'],unit=(u.deg,u.deg))
    raString = coor.ra.to_string(u.hourangle,precision=0)
    decString = coor.dec.to_string(u.deg,precision=0)
    ## Shorter way to name them is from g-r color
    gName, magString = [], []
    for ind, oneRow in enumerate(dat):
        gName.append("G-{:02d}-g-r={:.2f}".format(ind+1,oneRow['G_PS']-oneRow['R_PS']))
        #magString.append("{:.3f}".format(oneRow['G_PS']))
    t['Name'] =  gName
    tAlign['Name'] = alignName
    
    ## Assign priority from color
    t['Priority'] = 1000
    grColor = dat['G_PS'] - dat['R_PS']
    topPriority = np.abs(grColor - grSolarColor) < 0.02
    t['Priority'][topPriority] = 5000
    tAlign['Priority'] = -2
    
    ## Put in the Pan-Starrs G magnitudes
    t['Mag'] = dat['G_PS']
    tAlign['Mag'] = alignDat['G_PS']
    
    ## Put in the Gaia coordinates
    t['Coord Gaia'] = coor.to_string('hmsdms',sep=' ',precision=4)
    
    coorAlign = SkyCoord(alignDat['RA_ICRS'],alignDat['DE_ICRS'],unit=(u.deg,u.deg))
    tAlign['Coord Gaia'] = coorAlign.to_string('hmsdms',sep=' ',precision=4)
    
    ## Put in the Gaia epoch
    t['Epoch'] = 2015.5
    tAlign['Epoch'] = 2015.5
    
    ## Put in the Gaia equinox
    t['Equinox'] = 2000.0
    tAlign['Equinox'] = 2000.0
    
    ## Put in proper motion, it takes arcsec/yr
    t['pmRA'] = np.round(dat['pmRA'] / 1000.,4)
    t['pmDE'] = np.round(dat['pmDE'] / 1000.,4)
    tAlign['pmRA'] = np.round(alignDat['pmRA'] / 1000.,4)
    tAlign['pmDE'] = np.round(alignDat['pmDE'] / 1000.,4)
    #print(t)
    
    ## If second round, eliminate all sources already in mask
    if runRound == 1:
        roundText = ""
    if runRound == 2:
        maskColumns = ["X_STAR","Y_STAR","MIN_Y","MAX_Y","Percent","CCD_X","CCD_Y","Priority","Name","Mag"]
        prevMask = ascii.read('lists/autoslit/ngc2506_out.mask',names=maskColumns,
                             format='fixed_width_no_header',delimiter=' ')
        
        for oneRow in prevMask:
            matches = (t["Name"] == oneRow['Name'])
            t.remove_rows(matches)
        roundText = "_round2"
    
    t = astropy.table.vstack([tAlign,t])
    
    ## Find the center of the points
    coordCen = SkyCoord(np.mean(coor.ra),np.mean(coor.dec))
    coordCenString = coordCen.to_string('hmsdms',sep=' ',precision=4)
    
    t.insert_row(0,vals=["CENTER",9999,0.0,coordCenString,2015.0,2000.0,0.0,0.0])
    
    slitName = "ngc2506_solar_analogs{}.autoslit".format(roundText)
    
    t.write('lists/autoslit/{}'.format(slitName),
            format="ascii.fixed_width_no_header",
            delimiter=' ',overwrite=True)
    
    ## Also make the parameter file
    with open('lists/autoslit/input_par/ngc2506_solar_template.par') as templatePar:
        outPar = templatePar.readlines()#.append(templatePar.readline())
    
    outPar.append("BOXES {}\n".format(len(tAlign)))
    for oneBox in tAlign:
        outPar.append(str(oneBox['Name'])+"\n")
        
    outPar[0] = "FILENAME  {}\n".format(slitName)
    outPar[1] = "FILEOUT  ngc2506_out{}.mask\n".format(roundText)
    outPar.append('END\n')
    outPar.append('\n')
    
    with open('lists/autoslit/ngc2506_solar_analogs{}.par'.format(roundText),'w') as outParFile:
        outParFile.writelines(outPar)
    
    return t
    