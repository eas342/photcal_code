from astropy.io import fits,ascii
from astropy.coordinates import SkyCoord

dat = ascii.read('../lists_and_regions/coordinates_in_hms_dms_form.txt')
coor = SkyCoord(dat['Coordinates'])
outDat = dat
outDat['Ra (deg)'] = coor.ra.deg
outDat['Dec (deg)'] = coor.dec.deg
outDat.write('coordinates_in_hms_dms_and_degrees.csv',
             overwrite=True)
