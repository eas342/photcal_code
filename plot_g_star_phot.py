import matplotlib.pyplot as plt
from astropy.io import ascii

## Solar photometry, from http://mips.as.arizona.edu/~cnaw/sun.html
solDict = {}
solDict['g - r'] = 0.43
solDict['J mag - H mag'] = 3.65 - 3.33
solDict['K mag - IRAC1 mag'] = 3.27 - 3.26

gDat = ascii.read('g_star_photom.csv')

def plot_colmag(b1='g',b2='r',b3='J mag',b4='H mag'):
    """ Plot a color-color plot """
    
    fig, ax = plt.subplots()
    
    for oneType in ['G1','G2','G3']:
        tempPt = gDat['T Class'] == oneType
        ax.plot(gDat[tempPt][b1] - gDat[tempPt][b2],
                gDat[tempPt][b3] - gDat[tempPt][b4],'o',label=oneType)
        colorX = "{} - {}".format(b1,b2)
        ax.set_xlabel(colorX)
        colorY = "{} - {}".format(b3,b4)
        ax.set_ylabel(colorY)
    ax.axvline(solDict[colorX])
    ax.axhline(solDict[colorY])
    
    ax.legend()
    

def plot_colcol(b1='g',b2='r',b3='J mag',b4='H mag'):
    """ Plot a color-color plot """
    
    fig, ax = plt.subplots()
    
    for oneType in ['G1','G2','G3']:
        tempPt = gDat['T Class'] == oneType
        ax.plot(gDat[tempPt][b1] - gDat[tempPt][b2],
                gDat[tempPt][b3] - gDat[tempPt][b4],'o',label=oneType)
        colorX = "{} - {}".format(b1,b2)
        ax.set_xlabel(colorX)
        colorY = "{} - {}".format(b3,b4)
        ax.set_ylabel(colorY)
    ax.axvline(solDict[colorX])
    ax.axhline(solDict[colorY])
    
    ax.legend()
    
