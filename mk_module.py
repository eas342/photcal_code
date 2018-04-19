from astropy.io import ascii
from astropy.table import Table
import numpy as np
import hecto_module as hm
import phot_module as ps
import pdb
import os

class clusterClassification(object):
    """ Designed to get the data from the mkclass output """
    def __init__(self,mkOutFile='../classification/ngc2420_03_output.txt',src='NGC 2420'):
        self.mkOutFile = mkOutFile
        self.classData = ascii.read(mkOutFile,names=['SpFile','SpType','SpQuality','junk'])
        self.classData.sort('SpFile')
        self.src = src
        self.split_Types()
        className = os.path.basename(mkOutFile)
        ## The path to the synthesized photometry and classification
        self.synthPath = 'output_synthesized/class_phot_synth_{}.csv'.format(className)
    
    def split_Types(self):
        tclass, lumclass, xtra, tletter = [], [], [], []
        for oneType in self.classData['SpType']:
            splitArr = oneType.split(" ")
            numWords = len(splitArr)
            tclass.append(splitArr[0])
            tletter.append(splitArr[0][0])
            if numWords >= 2:
                lumclass.append(splitArr[1])
            else:
                lumclass.append('')
            if numWords >= 3:
                combExtra = r" ".join(splitArr[2:])
                xtra.append(combExtra)
            else:
                xtra.append(r"")
        self.classData['T Class'] = tclass
        self.classData['T Letter'] = tletter
        self.classData['Lum Class'] = lumclass
        self.classData['Extra Class info'] = xtra
    
    def save_sorted(self):
        baseNamePrefix = os.path.splitext(os.path.basename(self.mkOutFile))[0]
        self.classData.write('../classification/sorted/'+baseNamePrefix+'.csv')
    
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
        
    def get_phot_prev(self):
        """ Get previously saved photometry
        If file not found, it will generate one
        """
        if os.path.exists(self.synthPath) == False:
            self.get_phot(sType='All')
        self.photdat = ascii.read(self.synthPath)
    
    def get_phot(self,color='g-r',sType='G2 V',sCategory='SpType'):
        """ Gets the photometry of the stars of interest 
        Parameters
        -----------
        color: str
            The color to show (for example 'g-r')
        sType: str
            The spectral type to investigate. For example 'G2 V'
        sCategory: str
            The spectral category to search.
            For example, set sCategory = 'T Letter'and sTyp = 'G'
        """
        if self.src == 'NGC 2420':
            defaultIndices = [0,1]
        elif self.src == 'NGC 6811':
            defaultIndices = [0,1]
        else:
            defaultIndices = [0]
        
        hS = hm.clusterSpec(src=self.src,indices=defaultIndices)
        pS = ps.clusterPhot(src=self.src)
        uK = ps.clusterPhot(src=self.src,photType='UKIRTData')
        
        if (sType == None) | (sType == 'All'):
            lookRows = np.ones(len(self.classData),dtype=bool)
        else:
            lookRows = self.classData[sCategory] == sType
        t = Table()
        names, colors, mags, spTypeList = [], [], [], []
        TClassList, LClassList = [], []
        TLetterlist, xtraClassList = [], []
        posRA, posDec = [], []
        
        kmag, kmag_e = [], []
        
        for oneRow in self.classData[lookRows]:
            baseName = os.path.basename(oneRow['SpFile'])
            namePrefix = os.path.splitext(baseName)[0]
            searchName = namePrefix.split("_spec")[0]
            fibinfo = hS.getbyObjName(searchName)
            ra = float(fibinfo['RA'])
            dec = float(fibinfo['DEC'])
            posRA.append(ra)
            posDec.append(dec)
            phot = pS.lookup_src(ra,dec)
            
            names.append(str(fibinfo['OBJTYPE']))
            if phot is None:
                colors.append(np.nan)
                mags.append(np.nan)
            else:
                if ('g' in phot.colnames) & ('r' in phot.colnames):
                    colors.append(phot['g'] - phot['r'])
                    mags.append(phot['g'])
                else:
                    colors.append(np.nan)
                    mags.append(np.nan)
            
            uPhot = uK.lookup_src(ra,dec)
            if uPhot is None:
                kmag.append(np.nan)
                kmag_e.append(np.nan)
            else:
                kmag.append(uPhot['K mag'])
                kmag_e.append(uPhot['K mag err'])
            
            spTypeList.append(oneRow['SpType'])
            TClassList.append(oneRow['T Class'])
            LClassList.append(oneRow['Lum Class'])
            TLetterlist.append(oneRow['T Letter'])
            xtraClassList.append(oneRow['Extra Class info'])
            
        t['Name'] = names
        t['Color (g-r)'] = colors
        t['g'] = mags
        t['SpType'] = spTypeList
        t['T Class'] = TClassList
        t['Lum Class'] = LClassList
        t['T Letter'] = TLetterlist
        t['Extra Class info'] = xtraClassList
        t['ra'] = posRA
        t['dec'] = posDec
        t['K mag'] = kmag
        t['K mag err'] = kmag_e
        self.photdat = t
        
        if sType == 'All':
            t.write(self.synthPath,overwrite=True)

def get_candidates(t,widerT=True,gCut=1.0):
    """ 
    Get the interesting near-solar candidates from a 
    synthesized photometry + spectral classification table
    
    Parameters
    ----------
    t: Astropy Table
        A table of synthesized photometry + spectral classification
    widerT: bool
        Look at a wider range of temperature classes?
    gCut: float
        A magnitude cut in the g band.
        Eliminates all sources 1.0 magnitude brighter than the median
    """
    if widerT == True:
        tClasses = ['G0','G1','G2','G3']
    else:
        tClasses = ['G1','G2','G3']
    lClasses = ['IV','IV-V','V']
    
    TGood = np.zeros(len(t),dtype=np.bool)
    LGood = np.zeros(len(t),dtype=np.bool)
    
    for oneTClass in tClasses:
        goodP = (t['T Class'] == oneTClass)
        TGood = TGood | goodP
        
    for oneLClass in lClasses:
        goodP = (t['Lum Class'] == oneLClass)
        LGood = LGood | goodP
    
    medMag = np.median(t['g'][LGood & TGood])
    magCheck = t['g'] > (medMag - gCut)
    
    return t[LGood & TGood & magCheck]
    