# cluster_scaling_relations
Material related to Grudic &amp; Fall 2023 cluster scaling relations paper, including data compilation.

The compilation includes star cluster masses, ages, radii, and - where possible - environmental information such as galactocentric radii, Vcirc, tidal field, local gas, SFR, and stellar surface density, and gas velocity dispersion. Each cluster data row includes a `Reference` field containing the ADS bibcode for the works deserving credit for that dataset - *please cite them if you use that data point*. Some data unavailable in machine-readable format have been digitized from plots using WebPlotDigitizer; in these cases the data should be considered accurate to within a few per cent.

# Usage

`
from astropy.table import Table
Table.read("Cluster_Compilation.dat",format='ascii.basic')
`
