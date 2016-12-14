from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
import numpy as np

class clusterSpec(object):
    def __init__(self,
                redFile='../hectospec_data/2016.0207/reduction/0102/spHect-ngc2420_cat_1.2468-0102.fits'):
        hdulist = fits.open(redFile)
        self.fibinfo = Table(hdulist[5].data)
        self.makeCleanNames()
    
    def getbyObjName(self,objName):
        """ Search by object name
        """
        foundPts = self.fibinfo['cleanName'] == objName
        nFound = np.sum(foundPts)
        if np.sum(foundPts) == 1:
            return self.fibinfo[foundPts]
        else:
            print("Found "+str(nFound)+" objects for "+objName)
    
    def makeCleanNames(self):
        """
        Makes clean versions of the objtypes (no plus signs for use in file names)
        """
        cleanNames = []
        problemStrings = ['O-','-','+']
        replaceStrings = ['O_','_m','_p']
        for oneRow in self.fibinfo:
            newString = oneRow['OBJTYPE']
            for problemString, replaceString in zip(problemStrings,replaceStrings):
                newString = replaceString.join(newString.split(problemString))
            cleanNames.append(newString)
        self.fibinfo['cleanName'] = cleanNames