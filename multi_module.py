from astropy.io import ascii
from astropy.table import Table, vstack
from astropy.io import ascii
import numpy as np
import hecto_module as hm
import panstarrs_module as ps
import mk_module as mk
import pdb
import os
import pandas as pd
import yaml

def make_cluster_csv():
    """ Makes a CSV from the cluster file """
    dat = pd.read_excel('data/cluster_data.xlsx')
    dat.to_csv('data/cluster_data.csv',index=False)

def getRedFile(src,index=0,dataType='hectoData'):
    """ Finds the reduced file from the name and index """
    clusterF = yaml.load(open('data/cluster_files.yaml'))
    return clusterF[src][dataType][index]

def getClusterInfo(src,cProperty):
    """ Gets info on the cluster
    """
    cDat = ascii.read('data/cluster_data.csv')
    rowLook = cDat['Name'] == src
    return cDat[cProperty][rowLook]
    

def do_cm(fov=False,mkOutFile=None,src='NGC 2420',returnAx=False,
          figsize=None):
    """
    Does a color magnitude diagram
    """
    if mkOutFile is None:
        mkOutFile = getRedFile(src,dataType='mkClassification')
    
    psObj = ps.clusterPhot(src=src)
    if fov == True:
        psObj.plot_fov()
    else:
        psObj.plot_cm(figsize=figsize)
    
    colorShow = getClusterInfo(src,'g-r_solar')
    
    mkObj = mk.clusterClassification(mkOutFile=mkOutFile,src=src)
    typeExplore = ['F9 V','G0 V','G2 IV-V','G2 V','G4 V','G5 V','G8 V']
    colorArr = ['magenta','red','orange','green','olive','yellow','maroon']
    for oneType,dispCol in zip(typeExplore,colorArr):
        mkObj.get_phot(sType=oneType)
        if fov == True:
            psObj.ax.plot(mkObj.photdat['ra'],mkObj.photdat['dec'],'o',color=dispCol,
                          label=oneType)
        else:
            psObj.ax.plot(mkObj.photdat['Color (g-r)'],mkObj.photdat['g'],'o',color=dispCol,
                          label=oneType)
    if fov == False:
        psObj.ax.axvline(x=colorShow,linewidth=7.,alpha=0.3,color='red')
        
    psObj.ax.legend(loc='best',frameon=True)
    
    srcCleanName = src.replace(r" ",r"_")
    if fov == True:
        psObj.fig.savefig('plots/fov'+srcCleanName+'.pdf')
    else:
        custX=[-0.2,1.4]
        custY=[23,14]
        psObj.ax.set_xlim(custX[0],custX[1])
        psObj.ax.set_ylim(custY[0],custY[1])
        psObj.fig.savefig('plots/colormag'+srcCleanName+'.pdf')
    if returnAx == True:
        return psObj.fig, psObj.ax
    #psObj.fig.show()

def make_g2v_lists(sTypes=['G0 V','G2 V','G5 V']):
    """ Make a list of all near-G2V sources for all clusters"""
    clustFiles = yaml.load(open('data/cluster_files.yaml'))
    for oneClust in clustFiles.keys():
        mkOutFile = clustFiles[oneClust]['mkClassification']
        mkObj = mk.clusterClassification(mkOutFile=mkOutFile[0],src=oneClust)
        photDat = None
        for oneType in sTypes:
            #pdb.set_trace()
            mkObj.get_phot(sType=oneType)
            if len(mkObj.photdat) > 0:
                if photDat == None:
                    photDat = mkObj.photdat
                else:
                    photDat = vstack([photDat,mkObj.photdat])
        photDat.write('lists/g2vs'+oneClust+'.csv',overwrite=True)

def azProposal_plots():
    """
         Custom plots for Steward observatory proposal
    """
    fig, ax = do_cm(returnAx=True,figsize=(4,3.7))
    ax.set_xlim(0.1,0.9)
    ax.set_ylim(20,14)
    fig.savefig('plots/colormag_proposal.pdf',bbox_inches='tight')

def do_fov():
    """
    Shows the location of stars in the field 
    """
    do_cm(fov=True)
    