from astropy.io import fits,ascii
import pdb

dat = ascii.read('../lists/g_star_subset_ngc2506_ukirt.csv')

for oneColumn in dat.colnames:
    newN = oneColumn.replace(" ","_").replace("(","_").replace(")","_").replace("-","_")
    dat.rename_column(oneColumn,newN)

dat.write('../lists/g_star_subset_ngc2506_ukirt_ipac_form.tbl',
          format='ipac')
