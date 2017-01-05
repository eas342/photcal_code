from astropy.io import ascii
from astropy.table import Table
import numpy as np
import hecto_module as hm
import panstarrs_module as ps
import mk_module as mk
import pdb
import os

def do_cm():
    """
    Does a color magnitude diagram
    """
    psObj = ps.clusterPhot()
    psObj.plot_cm()
    mkObj = mk.clusterClassification()
    mkObj.get_phot()
    psObj.ax.plot(mkObj.photdat['Color (g-r)'],mkObj.photdat['g'],'o',color='red')
    psObj.fig.savefig('plots/colormag.pdf')
    #psObj.fig.show()