from astroquery.vizier import Vizier
from astropy.io import ascii
from astropy.table import Table, vstack, hstack
from astropy.io import ascii
import numpy as np
import pdb
import os
import pandas as pd
import yaml
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def plot_density_scatter(ax,all_x,all_y):
    
    # Calculate the point density
    goodpt = np.isfinite(all_x) & np.isfinite(all_y)
    x = all_x[goodpt]
    y = all_y[goodpt]
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    ax.scatter(x, y, c=z**4, s=10)
    ax.axhline(0.0,color='black',linestyle='dashed')

def get_2mass_cat(tab):
    """
    Get 2MASS catalog
    
    Parameters
    ----------
    tab: astropy table
        Table of objects to search. Must have "ra" and "dec" columns for sources
    """
    
    vizier = Vizier(catalog='II/246/out')
    vizier.ROW_LIMIT = 99999
    
    ra_med = np.median(tab['ra'])
    dec_med = np.median(tab['dec'])
    med_coord = SkyCoord(ra_med * u.deg, dec_med * u.deg)
    coord_tab = SkyCoord(tab['ra'] * u.deg,tab['dec'] * u.deg)
    sep = med_coord.separation(coord_tab)
    
    search_rad = np.max(sep) + 10. * u.arcsec
    
    tm_res = vizier.query_region(med_coord,radius=search_rad)[0]
    
    tm_coord = SkyCoord(tm_res['RAJ2000'],tm_res['DEJ2000'])
    idx, d2d, d3d = coord_tab.match_to_catalog_sky(tm_coord)
    
    max_sep = 1.0 * u.arcsec
    
    sep_constraint = d2d < max_sep
    bad_sep = d2d > max_sep
    
    for bandInd,oneBand in enumerate(['J','H','K']):
        
        tab['{}mag_2mass'.format(oneBand)] = tm_res['{}mag'.format(oneBand)][idx]
        tab['{}mag_2mass'.format(oneBand)][bad_sep] = np.nan
        
    return tab


def plot_twomass_ukirt(tab,src='NGC 2506'):
    fig, axArr = plt.subplots(3,2,figsize=(10,10))#,sharey=True)
    
    #axArr[0].plot(tab['J mag'],diff,'.')
    for bandInd,oneBand in enumerate(['J','H','K']):
    
        diff = tab['{} mag'.format(oneBand)] - tab['{}mag_2mass'.format(oneBand)]
        
        goodpt_diff = np.abs(diff) < 0.4 ## exlude outliers and bright points
        #goodpt_brightness = tab['K mag'] > 12.
        goodpt_brightness = (tab['K mag'] > 11.5) & (tab['K mag'] < 14.5)
        goodpt = goodpt_diff & goodpt_brightness
        medianOffset = np.nanmedian(diff[goodpt])
        print("Median UKIRT {} - 2MASS {} = {}".format(oneBand,oneBand,medianOffset))
        print("from  {} points".format(np.sum(goodpt)))
        axArr[bandInd,0].axhline(medianOffset,color='green',linestyle='dashed')
        axArr[bandInd,1].axhline(medianOffset,color='green',linestyle='dashed')
        
        densDat = plot_density_scatter(axArr[bandInd,0],tab['{} mag'.format(oneBand)],diff)
        axArr[bandInd,0].set_ylim(-0.2,0.2)
        axArr[bandInd,0].set_ylabel(r"{} (UKIRT) - {} (2MASS)".format(oneBand,oneBand))
        axArr[bandInd,0].set_xlabel("{} (UKIRT)".format(oneBand))
    
        densDat2 = plot_density_scatter(axArr[bandInd,1],tab['J mag'] - tab['K mag'],diff)
        axArr[bandInd,1].set_xlabel("J (UKIRT) - K (UKIRT)")
        axArr[bandInd,1].set_ylim(-0.2,0.2)
        axArr[bandInd,1].set_xlim(0,1)
    
    cax = fig.add_axes([0.91, 0.05, 0.02, 0.9]) #this locates the axis that is used for your colorbar. It is scaled 0 - 1.
    fig.colorbar(densDat2,cax,label='Point density')
#    minmax_show = [np.nanmin(tab['J mag']),np.nanmax(tab['Jmag_2mass'])]
    #plt.plot(minmax_show,1.0)
    outName = 'ukirt_2mass_{}.pdf'.format(src.replace(" ","_"))
    outPath = os.path.join('plots/ukirt_2mass_compare',outName)
    print("Saving plot to {}".format(outPath))
    fig.savefig(outPath,bbox_inches='tight')
    #fig.show()
"""
Vizier.query_object("NGC 2506",catalog="II/246/out")
result = Vizier.query_object("NGC 2506")

#catalog_list = Vizier.find_catalogs('NGC 2506')
tm_res = Vizier(catalog='II/246/out').query_region("NGC 2506",radius=0.3 *u.deg)[0]
catalog="II/246/out"
"""

