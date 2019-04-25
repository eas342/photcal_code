# coding: utf-8
imporfrom astropy.io import fits,ascii
from astropy.coordinates import SkyCoord
from astropy.io import fits,ascii
from astropy.coordinates import SkyCoord
dat = ascii.read('../lists_and_regions/coordinates_in_hms_dms_form.txt')
dat
dat = ascii.read('../lists_and_regions/coordinates_in_hms_dms_form.txt')
dat = ascii.read('../lists_and_regions/coordinates_in_hms_dms_form.txt')
dat
coor = SkyCoord(dat['Coordinates'])
coor
from copy import deepcopy
outDat = dat
outDat['Ra (deg)'] = coor.ra.deg
outDat['Dec (deg)'] = coor.dec.deg
outDat.write('coordinates_in_hms_dms_and_degrees.csv')
