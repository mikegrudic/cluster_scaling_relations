# cluster_scaling_relations
Material related to Grudic &amp; Fall 2023 cluster scaling relations paper, including data compilation.

The compilation includes star cluster masses, ages, radii, and - where possible - environmental information such as galactocentric radii, Vcirc, tidal field, local gas, SFR, and stellar surface density, and gas velocity dispersion. Each cluster data row includes a `Reference` field containing the ADS bibcode for the works deserving credit for that dataset - *please cite them if you use that data point*. Some data unavailable in machine-readable format have been digitized from plots using WebPlotDigitizer; in these cases the data should be considered accurate to within a few per cent.

# Usage

```
    from astropy.table import Table
    data = Table.read("Cluster_Compilation.dat",format='ascii.basic')
    M = data["Mass"]
    Reff = data["Reff"]
```

# Data fields

`Mass`: Cluster mass in solar
`Reff`: Cluster half-light radius in pc
`error_logReff`: Statistical error in half-light radius in dex (symmetrized if the original was asymmetric, as in Brown & Gnedin 2021)
`Rgc`: Galactocentric radius in kpc
`Region`: Galactic region label
`Vcirc`: Circular velocity at cluster galactocentric radius in km/s
`Omega`: Circular orbital frequency Vcirc/Rgc in km/s/kpc
