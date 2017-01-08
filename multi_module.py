from astropy.io import ascii
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import hecto_module as hm
import panstarrs_module as ps
import mk_module as mk
import pdb
import os
import pandas as pd

def make_cluster_file():
    """ Makes a CSV from the cluster file """
    dat = pd.read_excel('data/cluster_data.xlsx')
    dat.to_csv('data/cluster_data.csv',index=False)

def do_cm():
    """
    Does a color magnitude diagram
    """
    psObj = ps.clusterPhot()
    psObj.plot_cm()
    thisCluster = 'NGC 2420'
    cDat = ascii.read('data/cluster_data.csv')
    rowLook = cDat['Name'] == thisCluster
    colorShow = cDat['g-r_solar'][rowLook]
    
    mkObj = mk.clusterClassification()
    typeExplore = ['F9 V','G2 V','G5 V']
    colorArr = ['red','orange','green']
    for oneType,dispCol in zip(typeExplore,colorArr):
        mkObj.get_phot(sType=oneType)
        psObj.ax.plot(mkObj.photdat['Color (g-r)'],mkObj.photdat['g'],'o',color=dispCol,
                      label=oneType)
    
    psObj.ax.axvline(x=colorShow,linewidth=7.,alpha=0.3,color='red')
    psObj.ax.legend(loc='best')
    
    psObj.fig.savefig('plots/colormag.pdf')
    #psObj.fig.show()