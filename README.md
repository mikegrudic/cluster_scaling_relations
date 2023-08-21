# cluster_scaling_relations
Material related to Grudić &amp; Fall 2023 cluster scaling relations paper, including data compilation.

The compilation data file is `Cluster_Compilation.dat`. The compilation includes star cluster masses, ages, radii, and - where possible - environmental information such as galactocentric radii, Vcirc, epicyclic frequency, local gas, SFR, and stellar surface density, and gas velocity dispersion. Each cluster data row includes a `Reference` field containing the ADS bibcode for the works deserving credit for that dataset - *please cite them if you use that data point*. Some data unavailable in machine-readable format have been digitized from plots using WebPlotDigitizer; in these cases the data should be considered accurate to within a few per cent.

The data compilation procedure can be reproduced in Cluster_Compilation.py

# Example

```
    from astropy.table import Table
    data = Table.read("Cluster_Compilation.dat",format='ascii.basic')
    M_solar = data["Mass"]
    Reff_pc = data["Reff"]
```

# Data fields

`Reference`: ADS bibcodes of the papers that the cluster and galactic data came from - **PLEASE DON'T JUST CITE US; CITE THE NICE FOLKS WHO WORKED HARD GETTING THE DATA**

`Mass`: Cluster mass in solar

`Reff`: Cluster half-light radius in pc

`error_logReff`: Statistical error in half-light radius in dex (symmetrized if the original was asymmetric, as in Brown & Gnedin 2021)

`Rgc`: Galactocentric radius in kpc

`Region`: Galactic region label

`Vcirc`: Circular velocity at cluster galactocentric radius in km/s

`Omega`: Circular orbital frequency Vcirc/Rgc in km/s/kpc

`Kappa`: Epicyclic frequency in km/s/kpc

`Age`: Cluster age in yr

`SigmaGas`: Gas surface density in msun pc^-2

`SigmaStar`: Stellar surface density in msun pc^-2

`SigmaSFR`: SFR surface density in msun yr^-1 kpc^-2

`GalaxyDistance`: Distance to galaxy in kpc

`SigmaGas1D`: Plane-vertical gas velocity dispersion in km/s
