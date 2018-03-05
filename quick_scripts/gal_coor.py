# coding: utf-8
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii

dat = ascii.read('../data/cluster_data.csv')
coor = SkyCoord(dat['RA'] * u.degree,dat['Dec'] * u.degree)
#coor.galactic
print(coor.galactic.l)
print(coor.galactic.b)

