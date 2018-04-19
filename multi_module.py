from astropy.io import ascii
from astropy.table import Table, vstack
from astropy.io import ascii
import numpy as np
import hecto_module as hm
import phot_module as ps
import mk_module as mk
import pdb
import os
import pandas as pd
import yaml
from astropy.coordinates import SkyCoord
from astropy import units as u

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
    

def prep_all_spec():
    """ Prepares all Hectospec spectra for mkClass"""
    clustFiles = yaml.load(open('data/cluster_files.yaml'))
    cDat = ascii.read('data/cluster_data.csv')
    for cInd, oneClust in enumerate(clustFiles.keys()):
        nFiles = len(clustFiles[oneClust]['hectoData'])
        indSearch = np.arange(nFiles)
        hS = hm.clusterSpec(kernelWidth=cDat['Kernel Width'][cInd],
                            src=oneClust,indices=indSearch)
        hS.prepForMK()
    

def do_cm(fov=False,mkOutFile=None,src='NGC 2420',returnAx=False,
          figsize=None,photType='panStarrsData',
          groupType='near G2 V'):
    """
    Does a color magnitude diagram
    
    
    Parameters
    -----------
    fov: bool
          Show a field of view plot instead of the color magnitude diagram?
    mkOutFile: str
          Specify a different mkclass classification file.
          If None, then it will use te default for the object
    src: str
          The cluster source (e.g. "NGC 2506")
    returnAx: bool
          Return the axes object
    figsize: tuple
          The figure size as (x,y). If None, it will use defaults.
    photType: str
          The type of photometry 
    groupType: str
          How to group the stars.
          'All': show all stars w/ spectra
          'near G2 V': near G2V in spectral type
          'T Letter': F,G,K broadly
          'T Class': G0, G2, G5
          'Lum Class': V, IV-V, III, etc.
          'Candidates'
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
    
    if groupType == 'near G2 V':
        typeExplore = ['F9 V','G0 V','G2 IV-V','G2 V','G4 V','G5 V','G8 V']
        colorArr = ['magenta','red','orange','green','olive','yellow','maroon']
        sCategory = 'SpType'
    elif groupType in ['T Letter','T Class','Lum Class']:
        typeExplore = np.unique(mkObj.classData[groupType])
        colorArr = [None] * len(typeExplore)
        sCategory = groupType
    else:
        typeExplore = ['All']
        colorArr = ['red']
        sCategory = None
    
    mkObj.get_phot_prev()
    allPhotDat = mkObj.photdat
    
    if groupType == 'Candidates':
        if src == 'NGC 2420':
            widerT=False
        else:
            widerT=True
        ## Need a wider T range since we don't have as much data
        ## on NGC 2506 and NGC 6811 yet
        
        allPhotDat = mk.get_candidates(allPhotDat,widerT=widerT)
        typeExplore = [groupType]
    
    for oneType,dispCol in zip(typeExplore,colorArr):
        #mkObj.get_phot(sType=oneType,sCategory=sCategory)
        if sCategory is None:
            photDat = allPhotDat
        else:
            pts = allPhotDat[sCategory] == oneType
            photDat = allPhotDat[pts]
        if fov == True:
            psObj.ax.plot(photDat['ra'],photDat['dec'],'o',color=dispCol,
                          label=oneType)
        else:
            psObj.ax.plot(photDat['Color (g-r)'],photDat['g'],'o',color=dispCol,
                          label=oneType)
    if fov == False:
        psObj.ax.axvline(x=colorShow,linewidth=7.,alpha=0.3,color='red')
        
    psObj.ax.legend(loc='best',frameon=True)
    
    srcCleanName = src.replace(r" ",r"_")
    if fov == True:
        psObj.fig.savefig('plots/fov'+srcCleanName+'.pdf')
    else:
        custX=[-0.2,1.4]
        if src == 'NGC 6811':
            custY=[19,10]
        else:
            custY=[23,14]
        psObj.ax.set_xlim(custX[0],custX[1])
        psObj.ax.set_ylim(custY[0],custY[1])
        groupName = r"_".join(groupType.split(r" "))
        psObj.fig.savefig('plots/color_mag/colormag_{}_{}.pdf'.format(srcCleanName,groupName))
    if returnAx == True:
        return psObj.fig, psObj.ax
    #psObj.fig.show()

def cm_keck_proposal(src='NGC 2506'):
    """ Plots a color magnitude diagram for Keck proposal """
    srcCleanName = src.replace(r" ",r"_")
    
    p = ps.clusterPhot(src='NGC 2506')
    
    colorShow = getClusterInfo(src,'g-r_solar')
    p.plot_cm()
    p.ax.axvline(x=colorShow,linewidth=7.,alpha=0.3,color='red')
    
    coorFile = '../pan_starrs/pro/output/lris_targsNGC2506.csv'
    dat = ascii.read(coorFile)
    pts = dat['GROUP'] == 1
    p.ax.scatter(dat['G_M_R'][pts],dat['G'][pts],zorder=10,
                facecolors='none',edgecolors='red')
    
    p.ax.set_xlim(0,1.0)
    p.ax.set_ylim(21,14)
    
    p.fig.set_figwidth(5)
    p.fig.set_figheight(4)
    
    p.fig.savefig('plots/color_mag/colormag_selection_'+srcCleanName+'.pdf')

def make_g2v_lists(sTypes=['G0 V','G2 V','G5 V']):
    """ Make a list of all near-G2V sources for all clusters"""
    clustFiles = yaml.load(open('data/cluster_files.yaml'))
    for oneClust in clustFiles.keys():
        mkOutFile = clustFiles[oneClust]['mkClassification']
        mkObj = mk.clusterClassification(mkOutFile=mkOutFile[0],src=oneClust)
        mkObj.get_phot_prev()
        
        if oneClust == 'NGC 2420':
            widerT=False
        else:
            widerT=True
        
        photDat = mk.get_candidates(mkObj.photdat,widerT=widerT)
        
        coord = SkyCoord(ra=photDat['ra'] * u.degree,dec=photDat['dec'] * u.degree)
        photDat['coord hms'] = coord.to_string('hmsdms',sep=' ')
        photDat.write('lists/g_stars'+oneClust+'.csv',overwrite=True)

def azProposal_plots():
    """
         Custom plots for Steward observatory proposal
    """
    fig, ax = do_cm(returnAx=True,figsize=(4,3.7))
    ax.set_xlim(0.1,0.9)
    ax.set_ylim(20,14)
    fig.savefig('plots/for_proposals/colormag_proposal.pdf',bbox_inches='tight')

def do_fov():
    """
    Shows the location of stars in the field 
    """
    do_cm(fov=True)
    