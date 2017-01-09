from astropy.io import ascii
from astropy.table import Table
import numpy as np
import hecto_module as hm
import panstarrs_module as ps
import pdb
import os

class clusterClassification(object):
    """ Designed to get the data from the mkclass output """
    def __init__(self,mkOutFile='../classification/ngc2420_01_corrected_names.txt'):
        self.classData = ascii.read(mkOutFile,names=['SpFile','SpType','SpQuality','junk'])
        
    def get_counts(self):
        """
        Make a table of each spectral type and how many stars fit that spectral type
        Exampe: cMK = mk_module.clusterClassification() and cMK.get_counts()
        """
        uniq = np.unique(self.classData['SpType'])
        spTypes, counts = [], []
        
        for oneType in uniq:
            sameType = oneType == self.classData['SpType']
            spTypes.append(oneType)
            counts.append(np.sum(sameType))
        
        self.countTable = Table()
        self.countTable['SpType'] = spTypes
        self.countTable['Counts'] = counts
        
        self.countTable.pprint(max_lines=-1)
    
    def get_phot(self,color='g-r',sType='G2 V'):
        """ Gets the photometry of the stars of interest """
        hS = hm.clusterSpec()
        pS = ps.clusterPhot()
        if sType == None:
            lookRows = np.ones(len(self.classData),dtype=bool)
        else:
            lookRows = self.classData['SpType'] == sType
        t = Table()
        names, colors, mags, spTypeList = [], [], [], []
        posRA, posDec = [], []
        for oneRow in self.classData[lookRows]:
            baseName = os.path.basename(oneRow['SpFile'])
            namePrefix = os.path.splitext(baseName)[0]
            fibinfo = hS.getbyObjName(namePrefix)
            ra = fibinfo['RA']
            dec = fibinfo['DEC']
            posRA.append(ra)
            posDec.append(dec)
            phot = pS.lookup_src(ra,dec)
            names.append(fibinfo['OBJTYPE'])
            colors.append(phot['g'] - phot['r'])
            mags.append(phot['g'])
            spTypeList.append(oneRow['SpType'])
            
        t['Name'] = names
        t['Color (g-r)'] = colors
        t['g'] = mags
        t['SpType'] = spTypeList
        t['ra'] = posRA
        t['dec'] = posDec
        self.photdat = t
    