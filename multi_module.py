from astropy.io import ascii
from astropy.table import Table, vstack, hstack
from astropy.io import ascii
import numpy as np
import hecto_module as hm
import phot_module as ps
import mk_module as mk
import twomass_module as tm
import pdb
import os
import pandas as pd
import yaml
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt

def make_cluster_csv():
    """ Makes a CSV from the cluster file """
    dat = pd.read_excel('data/cluster_data.xlsx')
    dat.to_csv('data/cluster_data.csv',index=False)

def getRedFile(src,index=0,dataType='hectoData'):
    """ Finds the reduced file from the name and index """
    clusterF = yaml.safe_load(open('data/cluster_files.yaml'))
    return clusterF[src][dataType][index]

def getClusterInfo(src,cProperty):
    """ Gets info on the cluster
    """
    cDat = ascii.read('data/cluster_data.csv')
    rowLook = cDat['Name'] == src
    return cDat[cProperty][rowLook]
    

def prep_all_spec():
    """ Prepares all Hectospec spectra for mkClass"""
    clustFiles = yaml.safe_load(open('data/cluster_files.yaml'))
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
          The type of photometry (e.g. "panStarrsData" vs "UKIRTData")
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
    
    
    psObj = ps.clusterPhot(src=src,photType=photType)
    if fov == True:
        psObj.plot_fov()
    else:
        if photType == 'panStarrsData':
            color1='g'; color2='r' ; mag='g'
            solarColor = 'g-r_solar'
        else:
            color1='J mag' ; color2='H mag' ; mag='J mag'
            solarColor = 'J-H_solar'
        
        plotName = "{}-{}".format(color1.replace(" ","_"),color2.replace(" ","_"))
        
        psObj.plot_cm(figsize=figsize,
                      color1=color1,color2=color2,mag=mag)
        colorShow = getClusterInfo(src,solarColor)
    
    mkObj = mk.clusterClassification(mkOutFile=mkOutFile,src=src)
    
    if groupType == 'near G2 V':
        ## a list of categories to group sources into
        typeExplore = ['G0 IV to V','G1 IV to V','G2 IV to V','G3 IV to V']
        ## The colors of the categories
        colorArr = ['magenta','red','orange','green']
        ## The column name for the parameter that will be grouped
        sCategory = 'T Group'
    elif groupType in ['T Letter','T Class','Lum Class']:
        typeExplore = np.unique(mkObj.classData[groupType])
        colorArr = [None] * len(typeExplore)
        sCategory = groupType
    elif groupType == 'Candidates':
        typeExplore = ['Candidates']
        colorArr = ['red']
        sCategory = None
    else:
        typeExplore = ['All']
        colorArr = ['red']
        sCategory = None
    
    mkObj.get_phot_prev()
    allPhotDat = mkObj.photdat
    
    if (groupType == 'Candidates') | (groupType == 'near G2 V'):
        if src == 'NGC 2420':
            widerT=False
        else:
            widerT=True
        ## Need a wider T range since we don't have as much data
        ## on NGC 2506 and NGC 6811 yet
        
        allPhotDat = mk.get_candidates(allPhotDat,widerT=widerT)
    
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
            if photType == 'panStarrsData':
                xShow = photDat['Color (g-r)']
                yShow = photDat['g']
            else:
                xShow = photDat['J mag'] - photDat['H mag']
                yShow = photDat['J mag']
            
            psObj.ax.plot(xShow,yShow,'o',color=dispCol,
                          label=oneType)
    if fov == False:
        psObj.ax.axvline(x=colorShow,linewidth=7.,alpha=0.3,color='red')
        
    psObj.ax.legend(loc='best',frameon=True)
    psObj.ax.set_title(src)
    
    srcCleanName = src.replace(r" ",r"_")
    if fov == True:
        psObj.fig.savefig('plots/fov'+srcCleanName+'.pdf',bbox_inches='tight')
    else:
        if photType == 'panStarrsData':
            custX=[-0.2,1.4]
            if src == 'NGC 6811':
                custY=[19,10]
            else:
                custY=[23,14]
        else:
            custX = [0.1,0.8]
            custY = [19,11]
        psObj.ax.set_xlim(custX[0],custX[1])
        psObj.ax.set_ylim(custY[0],custY[1])
        groupName = r"_".join(groupType.split(r" "))
        psObj.fig.savefig('plots/color_mag/colormag_{}_{}_{}.pdf'.format(srcCleanName,groupName,plotName))
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
    clustFiles = yaml.safe_load(open('data/cluster_files.yaml'))
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
    
def autoslit_ukirt():
    """
    Find the UKIRT magnitudes of the LRIS slit mask stars put into autoslit
    
    """
    dat = Table.read('lists/gaia_coord/gaia_lris_targsNGC2506.fits')
    photDat =  ps.clusterPhot(src="NGC 2506",photType="UKIRTData")
    Jmag, Hmag, Kmag = [], [], []
    JmagErr, HmagErr, KmagErr = [], [], []
    for oneRow in dat:
        thisPhot = photDat.lookup_src(oneRow["RA_ICRS"],oneRow["DE_ICRS"])
        if thisPhot is None:
            Jmag.append(np.nan) ; JmagErr.append(np.nan)
            Hmag.append(np.nan) ; HmagErr.append(np.nan)
            Kmag.append(np.nan) ; KmagErr.append(np.nan)
        else:
            Jmag.append(thisPhot["J mag"]) ; JmagErr.append(thisPhot["J mag err"])
            Hmag.append(thisPhot["H mag"]) ; HmagErr.append(thisPhot["H mag err"])
            Kmag.append(thisPhot["K mag"]) ; KmagErr.append(thisPhot["K mag err"])
    
    dat['J mag'] = Jmag ; dat['J mag err'] = JmagErr
    dat['H mag'] = Hmag ; dat['H mag err'] = HmagErr
    dat['K mag'] = Kmag ; dat['K mag err'] = KmagErr
    dat['J - K mag'] = dat['J mag'] - dat['K mag']
    
    autoslitInfoBasename = 'lists/autoslit/candidate_info/ngc2506_autoslit'
    dat.write(autoslitInfoBasename+'.fits',overwrite=True)
    dat.write(autoslitInfoBasename+'.csv',overwrite=True)

def quick_photometry_ngc2506(specTypes=''):
    """
    Quick crude look at the stars from Keck LRIS
    """
    filePath = 'lists/sptype_002_solar_analogs/lists/early_lists/g_star_subset_ngc2506.csv'
    dat = ascii.read(filePath)
    phot = ps.clusterPhot(src='NGC 2506',photType='UKIRTData')
    jmag, hmag, kmag = [], [], []
    raDeg, decDeg = [], []
    for oneCoord in dat:
        coord = SkyCoord(oneCoord['RA'],oneCoord['Dec'],unit=(u.hourangle,u.deg))
        UKIRTphot = phot.lookup_src(ra=coord.ra.deg,dec=coord.dec.deg)
        if UKIRTphot == None:
            jmag.append('nan')
            hmag.append('nan')
            kmag.append('nan')
        else:
            jmag.append(UKIRTphot['J mag'])
            hmag.append(UKIRTphot['H mag'])
            kmag.append(UKIRTphot['K mag'])
        raDeg.append(coord.ra.deg)
        decDeg.append(coord.dec.deg)
    
    dat['J mag'] = jmag
    dat['H mag'] = hmag
    dat['K mag'] = kmag
    
    dat['RA (deg)'] = raDeg
    dat['Dec (deg)'] = decDeg
    
    outName = 'lists/sptype_002_solar_analogs/lists/g_star_subset_ngc2506_ukirt.csv'
    dat.write(outName,overwrite=True)
    

def make_coorFile():
    """ Makes a coordinate file for Rafia's Spitzer Pipeline
    """
    exampleCoorFile = '../spitzer/ngc2420_all_AOR_phot/ngc2420_example_skyCrd.csv'
    exCoorDat = ascii.read(exampleCoorFile)
    filePaths = np.unique(exCoorDat['AOR Filepath'])
    
    gDat = ascii.read("lists/g_starsNGC 2420.csv")
    coorDat = Table()
    
    paths, coors = [], []
    for oneFile in filePaths:
        for oneTarg in gDat:
            paths.append(oneFile)
            coorString = "{} {}".format(oneTarg['ra'],oneTarg['dec'])
            coors.append(coorString)
    coorDat['AOR Filepath'] = paths
    coorDat['Target Coordinates'] = coors
    
    print(coorDat)
    coorDat.write('lists/irac_coorFiles/ngc2420_sky2.csv')
    

def ukirtSpitzer():
    """
    Compare the Spitzer and IRAC photometry
    This is a check to see if the output from Rafia's pipeline looks reasonable
    """
    ## G star data with UKIRT and Pan Starrs info
    gDat = ascii.read("lists/g_starsNGC 2420.csv")
    coorG = SkyCoord(gDat['ra'],gDat['dec'],unit=(u.deg,u.deg))
    
    ## IRAC data
    #iDat = ascii.read('../spitzer/ngc2420_all_AOR_phot/ngc2420_final_aor_data.csv')
    iDat = ascii.read('../spitzer/ngc2420_all_AOR_phot/ngc2420_coor2_aor_data.csv')
    coorI = SkyCoord(iDat['Target RA'],iDat['Target Dec'],unit=(u.deg,u.deg))
    
    band = 'IRAC1'
    
    ## Get zero point and error
    if band == 'IRAC1':
        ## From here on 2019-01-23
        ## https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/17/
        zeroPoint = 280.9e3 ## mJy
        zeroErr = 4.1e3 ## mJy
    else:
        print("Unrecognized band!")
        return
    
    magZErr = zeroErr/zeroPoint * 2.5 / np.log(10.)
    
    iracMags, iracErrs = [], []
    iracFlux, iracFluxErr = [], []
    
    for indSrc, src in enumerate(gDat):
        indI = ((coorG[indSrc].separation(coorI) < 0.5 * u.arcsec) & 
                (iDat['FTime (sec)'] > 5.) & 
                (np.isfinite(iDat['Flux (mJy)'])))
        medFlux = np.median(iDat[indI]['Flux (mJy)'])
        standErr = np.median(iDat[indI]['Error (mJy)']) / np.sum(indI)
        
        medMag = -2.5 * np.log10(medFlux/zeroPoint)
        medMagErr = standErr / medFlux * 2.5  / np.log(10.)
        totErr = np.sqrt(medMagErr**2 + magZErr**2)
        
        iracMags.append(medMag)
        iracErrs.append(totErr)
        iracFlux.append(medFlux)
        iracFluxErr.append(standErr)
    
    gDat[band+ ' mag'] = iracMags
    gDat[band+' err'] = iracErrs
    gDat[band+' Flux (mJy)'] = iracFlux
    gDat[band+' Flux Err (mJy)'] = iracFluxErr
    
    ## Also throw in Pan Starrs Data
    photDat = ps.clusterPhot(src="NGC 2420",photType="panStarrsData")
    
    for ind, oneRow in enumerate(gDat):
        thisPhot = photDat.lookup_src(oneRow['ra'],oneRow['dec'])
        if ind == 0:
            allPSphot = thisPhot
        else:
            allPSphot = vstack([allPSphot,thisPhot])
    
    gDat = hstack([gDat,allPSphot[['dg','r','dr','i','di','z','dz']]])
    
    return gDat
    

def example_gaia_source_aptx_file():
    """
    Look at the color magnitude diagram for the example source
    """
    cPhot = ps.clusterPhot(src='NGC 2506',photType='UKIRTData')
    
    coord = SkyCoord('07 59 58.5995 -10 46 46.93',unit=(u.hourangle,u.deg))
    
    g1 = cPhot.lookup_src(ra=coord.ra.deg,dec=coord.dec.deg)
    
    cPhot.plot_cm(color1='J mag',color2='K mag',mag='J mag')
    cPhot.ax.plot(g1['J mag'] - g1['K mag'],g1['J mag'],'o',markersize=16,color='red')
    cPhot.ax.text(g1['J mag'] - g1['K mag'],g1['J mag'],'GAIA DR2 3038044777529190144',color='red')
    
    cPhot.ax.set_xlim(0.15,0.9)
    cPhot.ax.set_ylim(18,12)
    cPhot.fig.savefig('plots/specific_plots/gaia_dr_3038044777529190144.png')
    plt.close(cPhot.fig)

def cm_autoslit(color='J - K'):
    """
    Plot the color magnitude of the LRIS slit mask stars put into autoslit
    
    Parameters
    -------------
    color: str
        Specify which color to plot. E.g. "J - K" or "g - r"
    
    """
    autoDat = Table.read("lists/autoslit/candidate_info/ngc2506_autoslit.fits")
    if color == 'J - K':
        photColor = autoDat['J - K mag']
        colorName = 'J - K'
        ## using Christopher Willmer's Abs Mag of the sun page
        solarColor = 5.03 - 4.64
    elif color == 'J - H':
        photColor = autoDat['J mag'] - autoDat['H mag']
        colorName = 'J - H'
        ## using Christopher Willmer's Abs Mag of the sun page
        solarColor = 0.32
    elif color == 'g - r':
        photColor = autoDat['G_M_R_PS']
        colorName = 'g - r'
        ## using my notebook with extinction
        solarColor = 0.43
    else:
        print("Unrecognized color")
        return
    
    mag = autoDat['K mag']
    magName = 'K'
    
    fig, ax = plt.subplots()
    ax.plot(photColor,mag,'.')
    ax.invert_yaxis()
    
    ax.axvline(solarColor)
    ax.set_xlabel(colorName)
    ax.set_ylabel(magName)
    fig.savefig('plots/autoslit/{}.pdf'.format(color.replace(" ","_")))

def compare_ukirt_2mass(src='NGC 2506'):
    """
    Compare UKIRT and 2MASS
    """
    cPhot  = ps.clusterPhot(src=src,photType='UKIRTData')
    tm_res = tm.get_2mass_cat(cPhot.dat)
    tm.plot_twomass_ukirt(tm_res,src=src)
    