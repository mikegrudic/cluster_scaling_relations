# %%

from astropy.table import Table
from os import system
from os.path import exists
from astropy.coordinates import ICRS, Distance, Angle
from astropy.table import Table
from astropy import units as u
import numpy as np

include_galaxies_without_env_data = True


def correct_rgc(
    coord,
    glx_ctr=ICRS("00h42m44.33s", "+41d16m07.5s"),
    glx_PA=Angle("37d42m54s"),
    glx_incl=Angle("77.5d"),
    glx_dist=Distance(783, unit=u.kpc),
):
    """Computes deprojected galactocentric distance.
    Inspired by: http://idl-moustakas.googlecode.com/svn-history/
        r560/trunk/impro/hiiregions/im_hiiregion_deproject.pro
    Parameters
    ----------
    coord : :class:`astropy.coordinates.ICRS`
        Coordinate of points to compute galactocentric distance for.
        Can be either a single coordinate, or array of coordinates.
    glx_ctr : :class:`astropy.coordinates.ICRS`
        Galaxy center
    glx_PA : :class:`astropy.coordinates.Angle`
        Position angle of galaxy disk.
    glx_incl : :class:`astropy.coordinates.Angle`
        Inclination angle of the galaxy disk.
    glx_dist : :class:`astropy.coordinates.Distance`
        Distance to galaxy.
    Returns
    -------
    obj_dist : class:`astropy.coordinates.Distance`
        Galactocentric distance(s) for coordinate point(s).
    """
    # distance from coord to glx centre
    sky_radius = glx_ctr.separation(coord)
    avg_dec = 0.5 * (glx_ctr.dec + coord.dec).radian
    x = (glx_ctr.ra - coord.ra) * np.cos(avg_dec)
    y = glx_ctr.dec - coord.dec
    # azimuthal angle from coord to glx  -- not completely happy with this
    phi = glx_PA - Angle("90d") + Angle(np.arctan(y.arcsec / x.arcsec), unit=u.rad)

    # convert to coordinates in rotated frame, where y-axis is galaxy major
    # ax; have to convert to arcmin b/c can't do sqrt(x^2+y^2) when x and y
    # are angles
    xp = (sky_radius * np.cos(phi.radian)).arcmin
    yp = (sky_radius * np.sin(phi.radian)).arcmin

    # de-project
    ypp = yp / np.cos(glx_incl.radian)
    obj_radius = np.sqrt(xp**2 + ypp**2)  # in arcmin
    obj_dist = Distance(
        Angle(obj_radius, unit=u.arcmin).radian * glx_dist, unit=glx_dist.unit
    )

    # Computing PA in disk (unused)
    obj_phi = Angle(np.arctan(ypp / xp), unit=u.rad)
    # TODO Zero out very small angles, i.e.
    # if np.abs(Angle(xp, unit=u.arcmin)) < Angle(1e-5, unit=u.rad):
    #     obj_phi = Angle(0.0)

    return obj_dist


ids = []
galaxy_name = []  # galaxy name
references = []  # bibcodes for references the data came from
masses = []  # best mass estimate
# mass_upper = [] # +1sigma mass quantile if available
# mass_lower = [] # -1sigma mass quantile if available
reffs = []  # projected effective (half light) radius
reff_error = []  # symmetrized 1sigma RELATIVE error in dex
# reff_lower = [] # -1sigma effective radius quantile if available
Rgcs = []  # galactocentric radius
regions = []  # galactic subregion designation (if any)
rotationcurve_data = []  # data file containing rotation curve data
Vcs = []  # circular velocity at Rgc
omegas = []  # orbital angular frequency = Vc/Rgc
kappas = []  # epicyclic frequency
rho_tidals = (
    []
)  # tidal density = 3 M / (4 pi r_J) where r_J = (G M / (4 omega^2 - kappa^2))
sigma_gass = []  # gas surface density
sigma_SFRs = []  # SFR surface density
sigma_stars = []  # stellar surface density
DGCs = []  # distance to galaxy in Mpc
ages = []
sigma_1Ds = []

compilation = Table()

# %% [markdown]
# # M82 clusters from Mayya 2008; Cuevas 2019,2021

# %%
# from 2021MNRAS.500.4422C supplementary files
id = np.genfromtxt("tabla_1.dat", dtype=str, skip_header=1)[:, 0]
data = np.genfromtxt("tabla_1.dat", skip_header=1)
RGC = data[:, 1]
Rh = data[:, 8]
Rh_error = (data[:, 9] + data[:, 10]) / 2
Rh_error = Rh_error / Rh / np.log(10)
Rh_king = data[:,]
M = 10 ** data[:, 21]
# print(M,Rh)
# plt.loglog(M,Rh,ls="",marker='o')

# model of M82 following Schneider 2018, fitting rotation curve to Greco 2012
Rgc2, Vrot = np.loadtxt("M82_RGC_vs_Vrot.csv").T
Menc = Vrot**2 * Rgc2 * 2.32e5
omega = Vrot / Rgc2
kappa = np.sqrt(2 * omega / Rgc2 * np.gradient(Rgc2**2 * omega) / np.gradient(Rgc2))
rho_tidal = 3 * (4 * omega**2 - kappa**2) / (4 * np.pi) * 2.325e-4

from scipy.integrate import trapz

Rgc3 = np.logspace(-1, 1, 100)
R_stellar = 0.8
M_stellar = 6e9
M_gas = 4e9
R_gas = 2 * R_stellar
sigma_gas = np.exp(-Rgc3 / R_gas)
sigma_gas *= M_gas / trapz(2 * np.pi * Rgc3 * sigma_gas, Rgc3) / 1e6
# plt.plot(Rgc3, sigma_gas); plt.yscale('log')
sigma_gas = np.interp(RGC, Rgc3, sigma_gas)
sigma_star = np.exp(-Rgc3 / R_stellar)
sigma_star *= M_stellar / trapz(2 * np.pi * Rgc3 * sigma_star, Rgc3) / 1e6
sigma_star = np.interp(RGC, Rgc3, sigma_star)
sigma_SFR = np.repeat(np.nan, len(M))


Ncl = len(M)
galaxy_name.append(np.repeat("M82", Ncl))
ids.append(id)
references.append(
    np.repeat(
        "2008ApJ...679..404M;2020MNRAS.492..993C;2021MNRAS.500.4422C;2012ApJ...757...24G",
        Ncl,
    )
)  # bibcodes for references the data came from
masses.append(M)  # best mass estimate
reffs.append(Rh)  # projected effective (half light) radius
reff_error.append(Rh_error)
Rgcs.append(RGC)  # galactocentric radius
regions.append(np.repeat("", len(M)))  # galactic subregion designation (if any)
Vcs.append(np.interp(RGC, Rgc2, Vrot))  # circular velocity at Rgc
omegas.append(np.interp(RGC, Rgc2, omega))  # orbital angular frequency = Vc/Rgc
kappas.append(np.interp(RGC, Rgc2, kappa))  # epicyclic frequency
rho_tidals.append(np.interp(RGC, Rgc2, rho_tidal))
sigma_gass.append(sigma_gas)
sigma_SFRs.append(sigma_SFR)
sigma_stars.append(sigma_star)
ages.append(
    np.array([(1.00001e7 if r < 0.4 else 1e8) for r in RGC])
)  #  VERY approximate age, do not use for analysis - just saying stuff in the central region is young (~10Myr), while stuff in the disk is ~100Myr old
DGCs.append(np.repeat(3630, Ncl))  # Freedman 1994, as assumed in Cuevas-Otahola 2020;21
sigma_1Ds.append(
    np.repeat(np.nan, Ncl)
)  # gas velocity dispersion profile not available; inclination too high

# %% [markdown]
# # M31 PHAT clusters
# Here we cross-reference effective projected radii from 2015ApJ...802..127J and masses/ages from 2016ApJ...827...33J

# %%
system(
    "wget -N https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJ/802/127/table2.dat"
)
t1 = Table.read("fits?J%2FApJ%2F802%2F127%2Ftable2.dat")
system(
    "wget -N https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJ/827/33/table4.dat"
)
t2 = Table.read("fits?J%2FApJ%2F827%2F33%2Ftable4.dat")
dist = 761  # distance of 761+/-11 kpc - 2021ApJ...920...84L
id1 = t1["APID"]  # ID from catalog 1
id2 = t2["APID"]  # ID from catalog 2
mag = t1["F475W"]
mag = np.array([dict(zip(id1, mag))[i] for i in id2])
cut = mag < 100  # 20.5)
region = np.array([str(t).split()[0] for t in t2["RegID"]])  # region IDs from catalog 2
RA = t1["RAdeg"]  # right ascension from catalog 2
DEC = t1["DEdeg"]  # declination from catalog
coords = [
    ICRS("%gd" % RA[id1 == i], "%gd" % DEC[id1 == i]) for i in id2
]  # deprojected galactocentric radius
Rgc = np.array([correct_rgc(c).value for c in coords])
age = 10 ** t2["logAge"]
mass = 10 ** t2["logMass"]
m_min = 10 ** t2["b_logMass"]
m_max = 10 ** t2["B_logMass"]
reff = t1["Reff"] * dist * 1e3 * 4.484e-6
reff = np.array([dict(zip(id1, reff))[i] for i in id2])

print(np.median(reff))

# below are values reported in Johnson 2016 Tables 1 and 2 for the different regions
sigma_gasdict = {
    "1": 4.47,
    "2a": 10.02,
    "2": 10.54,
    "2b": 12.01,
    "2c": 11.34,
    "2d": 9.84,
    "2e": 9.66,
    "3": 5.34,
    "NA": np.nan,
}
sigma_SFRdict = {
    "1": 10**-2.96,
    "2a": 10**-2.45,
    "2": 10**-2.55,
    "2b": 10**-2.59,
    "2c": 10**-2.61,
    "2d": 10**-2.65,
    "2e": 10**-2.48,
    "3": 10**-3.13,
    "NA": np.nan,
}
sigma_1Ddict = {
    "1": 7.94,
    "2a": 9.73,
    "2": 8.65,
    "2b": 9.46,
    "2c": 8.03,
    "2d": 8.16,
    "2e": 8.12,
    "3": 7.17,
    "NA": np.nan,
}
# RGCdict = {'1': 6.61, '2a': 11.8, '2': 11.45, '2b': 12.14, '2c': 12.11, '2d': 12.16, '2e': 11.14, '3': 15.83, 'NA': np.nan}

sigma_gas = [sigma_gasdict[r] for r in region]  # for each cluster
sigma_SFR = [sigma_SFRdict[r] for r in region]  # for each cluster
sigma_1D = [sigma_1Ddict[r] for r in region]  # for each cluster
# STELLAR SURFACE DENSITY CALCULATION
# lift stellar bulge+disk mass model from Tamm 2012A&A...546A...4T Table A.1, these are fits to the ellipsoidal Einasto model, rho = rho_c exp(-d_N ((a/a_c)^(1/N)-1)), where a = sqrt(r^2 + z^2/q^2)
from scipy.integrate import quad

sersic_params_bulge = (2.025, 0.73, 4.0, 11.67, 2.2e-1 * 1e9, 4.9)
sersic_params_disk = (11.35, 0.1, 1.0, 2.67, 1.72e-2 * 1e9, 4.8)
rgrid = np.logspace(-1, 2, 1000)  # grid of radii
sigma_star = np.zeros_like(rgrid)


def sersic_profile(z, r, ac, q, N, d_N, rhoc):
    a = np.sqrt(r * r + z * z / (q * q))
    return rhoc * np.exp(-d_N * ((a / ac) ** (1.0 / N) - 1))


for p in (
    sersic_params_bulge,
    sersic_params_disk,
):  # compute surface density of each component and sum
    ac, q, N, d_N, rhoc, Mtot = p
    sigma_star += np.array(
        [
            quad(sersic_profile, -np.inf, np.inf, args=(R, ac, q, N, d_N, rhoc))[0]
            for R in rgrid
        ]
    )  # grid of radii, integrating along the z axis to get surface density

sigma_star /= 1e6  # convert to msun pc^-2

sigma_star = np.interp(
    Rgc, rgrid, sigma_star
)  # evaluate stellar surface density at cluster positions

Rgc2, Vrot = np.loadtxt("M31_R_vs_Vrot.csv").T

Menc = Vrot**2 * Rgc2 * 2.32e5
omega = Vrot / Rgc2
kappa = np.sqrt(2 * omega / Rgc2 * np.gradient(Rgc2**2 * omega) / np.gradient(Rgc2))
rho_tidal = (
    3 * (4 * omega**2 - kappa**2) / (4 * np.pi) * 2.325e-4
)  # 3 * Menc / (4*np.pi*(1e3*Rgc2)**3)
# plt.loglog(Rgc2,Vrot)

galaxy_name.append(np.repeat("M31", len(mass[cut])))
ids.append(id2[cut])
references.append(
    np.repeat(
        "2012A&A...546A...4T;2015ApJ...802..127J;2016ApJ...827...33J", len(mass[cut])
    )
)  # bibcodes for references the data came from
masses.append(mass[cut])  # best mass estimate
reffs.append(reff[cut])  # projected effective (half light) radius
Rgcs.append(Rgc[cut])  # galactocentric radius
regions.append(region[cut])  # galactic subregion designation (if any)
Vcs.append(np.interp(Rgc[cut], Rgc2, Vrot))  # circular velocity at Rgc
omegas.append(np.interp(Rgc[cut], Rgc2, omega))  # orbital angular frequency = Vc/Rgc
kappas.append(np.interp(Rgc[cut], Rgc2, kappa))  # epicyclic frequency
rho_tidals.append(np.interp(Rgc[cut], Rgc2, rho_tidal))
sigma_gass.append(np.array(sigma_gas)[cut])
sigma_SFRs.append(np.array(sigma_SFR)[cut])
sigma_stars.append(sigma_star[cut])
ages.append(age[cut])
DGCs.append(
    np.repeat(dist, cut.sum())
)  # distance of 761+/-11 kpc - 2021ApJ...920...84L
sigma_1Ds.append(np.array(sigma_1D)[cut])
reff_error.append(np.zeros_like(reffs[-1]))

# %% [markdown]
# # M33 PHATTER clusters
# Here we cross-reference effective projected radii from Johnson 2022 and masses/ages from Wainer 2022

# %%
t1 = Table.read("M33_Star_Cluster_Cat_with_CMD_Ests_2023Update_USE_THIS_ONE.fits")
t2 = Table.read("Final_M33_catalog_with_CMD_Estimates.fits")
dist = 875  # 2022ApJ...933..201L - 873 +/- 20kpc
id1 = t1["ID"]  # ID from catalog 1
id2 = t2["ID"]  # ID from catalog 2
mag = t1["MAG475"]
magdict = dict(zip(id1, mag))
mag = np.array([magdict[i] for i in id2])
cut = (
    mag < 100
)  # F475W magnitude cut if desired; supposedly incomplete below 20.5 - Johnson 2016
RA = t1["RA"]  # right ascension from catalog 2
DEC = t1["DEC"]  # declination from catalog
coords = [
    ICRS("%gd" % RA[id1 == i], "%gd" % DEC[id1 == i]) for i in id2
]  # deprojected galactocentric radius
Rgc = np.array(
    [
        correct_rgc(
            c,
            glx_ctr=ICRS("01h33m50.02s", "+30d39m36.7s"),
            glx_dist=Distance(870, unit=u.kpc),
            glx_PA=Angle("21d"),
            glx_incl=Angle("50d"),
        ).value
        for c in coords
    ]
)  # distance from Galleti 2004; de Grijs & Bono 2014. PA and inc from Corbelli 1999
# age = 10**t2['Age_Med']
# mass = 10**t2['Mass_Med']
age = 10 ** t2["Bestfit_Age"]
mass = 10 ** t2["Bestfit_Mass"]
reff = t1["RADIUS_EFF_ARCSEC"] * dist * 1e3 * 4.484e-6
reffdict = dict(zip(id1, reff))
reff = np.array([reffdict[i] for i in id2])

Rgc2, sigma_star = np.loadtxt("M33_RGC_vs_SigmaStar.csv").T

sigma_star = np.interp(
    Rgc, Rgc2, sigma_star
)  # evaluate stellar surface density at cluster positions

Rgc2, sigma_gas = np.loadtxt("M33_RGC_vs_SigmaGas.csv").T  # Gonzalez-Lopezlira 2012

sigma_gas = np.interp(
    Rgc, Rgc2, sigma_gas
)  # evaluate gas surface density at cluster positions

Rgc2, sigma_SFR = np.loadtxt("M33_RGC_vs_SigmaSFR.csv").T  # Gonzalez-Lopezlira 2012

sigma_SFR = np.interp(
    Rgc, Rgc2, sigma_SFR
)  # evaluate SFR surface density at cluster positions

Rgc2, sigma_1D = 2.909e-4 * dist * np.array(
    [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
), [
    8.5,
    8.2,
    8.3,
    7.6,
    7.7,
    6.8,
    7.4,
    7.5,
    6.7,
    6.6,
    6.4,
    8.0,
]  # Kam et al 2017AJ....154...41K, converting RGC from arcsec

sigma_1D = np.interp(
    Rgc, Rgc2, sigma_1D
)  # evaluate 1D gas velocity dispersion at cluster positions

Rgc2, Vrot = np.loadtxt("M33_RGC_vs_Vrot.csv").T

Menc = Vrot**2 * Rgc2 * 2.32e5
omega = Vrot / Rgc2
kappa = np.sqrt(2 * omega / Rgc2 * np.gradient(Rgc2**2 * omega) / np.gradient(Rgc2))
rho_tidal = (
    3 * (4 * omega**2 - kappa**2) / (4 * np.pi) * 2.325e-4
)  # 3 * Menc / (4*np.pi*(1e3*Rgc2)**3)

galaxy_name.append(np.repeat("M33", len(mass[cut])))
ids.append(id2[cut])
references.append(
    np.repeat(
        "2012ApJ...761..124G;2022ApJ...928...15W;2022ApJ...938...81J", len(mass[cut])
    )
)  # bibcodes for references the data came from
masses.append(mass[cut])  # best mass estimate
reffs.append(reff[cut])  # projected effective (half light) radius
Rgcs.append(Rgc[cut])  # galactocentric radius
regions.append(np.repeat("", cut.sum()))  # galactic subregion designation (if any)
Vcs.append(np.interp(Rgc[cut], Rgc2, Vrot))  # circular velocity at Rgc
omegas.append(np.interp(Rgc[cut], Rgc2, omega))  # orbital angular frequency = Vc/Rgc
kappas.append(np.interp(Rgc[cut], Rgc2, kappa))  # epicyclic frequency
rho_tidals.append(np.interp(Rgc[cut], Rgc2, rho_tidal))
sigma_gass.append(np.array(sigma_gas)[cut])
sigma_SFRs.append(np.array(sigma_SFR)[cut])
sigma_stars.append(sigma_star[cut])
ages.append(age[cut])
DGCs.append(np.repeat(dist, cut.sum()))  # 2022ApJ...933..201L - 875 +/- 20kpc
sigma_1Ds.append(np.array(sigma_1D)[cut])
reff_error.append(np.zeros_like(reffs[-1]))

# %% [markdown]
# # Large and Small Magellanic Clouds
# Data are from McLaughlin & van der Marel 2005; we adopt their King model fits. Rotation curve for SMC from di Teodoro 2018, LMC from Alves 2000

# %%
##### SMC and LMC ###########################################################################################################################
system(
    "wget -N https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJS/161/304/table11.dat"
)
system(
    "wget -N https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJS/161/304/table12.dat"
)
system(
    "wget -N https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJS/161/304/table14.dat"
)
system(
    "wget -N https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/fits?J/ApJS/161/304/table8.dat"
)

t1 = Table.read("J_ApJS_161_304_table11.dat.fits")
t2 = Table.read("fits?J%2FApJS%2F161%2F304%2Ftable8.dat")
t3 = Table.read("fits?J%2FApJS%2F161%2F304%2Ftable12.dat")
t4 = Table.read("fits?J%2FApJS%2F161%2F304%2Ftable14.dat")

galaxies = [n.split("-")[0] for n in t1["Cluster"]]

dist_dict = {"LMC": 49.59, "SMC": 62.8}  # Pieterynski 2019  # Cioni 2000

sigma_SFR_dict = {"LMC": 3.66e-3, "SMC": 1e-3}  # Baumgardt 2013  # Goddard 2010

r, menc = np.loadtxt("SMC_R_vs_Menc_stellar.csv").T  # Bekki 2009
sigma_star = np.gradient(menc) / np.gradient(np.pi * r**2) / 1e6
np.savetxt("SMC_R_vs_SigmaStar.csv", np.c_[r, sigma_star])

r, menc = np.loadtxt("SMC_R_vs_Menc_gas.csv").T  # Bekki 2009
sigma_gas = np.gradient(menc) / np.gradient(np.pi * r**2) / 1e6
np.savetxt("SMC_R_vs_SigmaGas.csv", np.c_[r, sigma_gas])

for g in np.unique(galaxies):
    if not exists(g + "_R_vs_Vrot.csv"):
        continue
    model = t1["Mod"]
    print(np.unique(model).data)
    mask = (
        (model == b"K ") * (np.array(galaxies) == g) * np.isfinite(t1["Rh"])
    )  # select King model fits
    ids.append(t1["Cluster"][mask])
    agedict = dict(zip(t2["Cluster"], 10 ** t2["logAge"]))
    ages.append(np.array([agedict[n] for n in t3["Cluster"][mask]]))
    mass = 10 ** t3["Mtot"][mask]
    masses.append(mass)
    reffs.append(10 ** t1["Rh"][mask])
    reff_error.append(0.5 * (t1["E_Rh"] + t1["e_Rh"])[mask])
    Rgcs.append(t4["Rad"][mask])
    Rgc = t4["Rad"][mask]
    Rgc2, Vrot = np.loadtxt(
        g + "_R_vs_Vrot.csv"
    ).T  # rotation curves from Alves 2000 for LMC and

    omega = Vrot / Rgc2
    kappa = np.sqrt(
        2 * omega / Rgc2 * np.gradient(Rgc2**2 * omega) / np.gradient(Rgc2)
    )
    rho_tidal = 3 * (4 * omega**2 - kappa**2) / (4 * np.pi) * 2.325e-4

    R_sigma_star, sigma_star = np.loadtxt(
        g + "_R_vs_SigmaStar.csv"
    ).T  # LMC data from Wong 2009; SMC data from Bekki 2009 NFW mass model
    sigma_star = np.interp(Rgc, R_sigma_star, sigma_star)

    if g == "LMC":
        R_sigma_gas, sigma_gas2 = np.loadtxt(
            "LMC_RGC_vs_SigmaGasAtomic.csv"
        ).T  # LMC data from from Wong 2009; SMC data from Bekki 2009 NFW mass model
        sigma_gas = np.interp(Rgc, R_sigma_gas, sigma_gas2)
        R_sigma_gas, sigma_gas2 = np.loadtxt(
            "LMC_RGC_vs_SigmaGasMolecular.csv"
        ).T  # LMC data from from Wong 2009; SMC data from Bekki 2009 NFW mass model
        sigma_gas += np.interp(Rgc, R_sigma_gas, sigma_gas2)
    else:
        R_sigma_gas, sigma_gas = np.loadtxt(
            g + "_R_vs_SigmaGas.csv"
        ).T  # SMC data from Bekki 2009 NFW mass model
        sigma_gas = np.interp(Rgc, R_sigma_gas, sigma_gas)

    if g == "SMC":
        R_sigma_1D, sigma_1D = np.loadtxt(
            g + "_R_vs_Sigma1D.csv"
        ).T  # SMC data from Di Teodoro 2019MNRAS.483..392D
        sigma_1D = np.interp(Rgc, R_sigma_1D, sigma_1D)

    Ncl = len(mass)
    galaxy_name.append(np.repeat(g, Ncl))
    references.append(
        np.repeat(
            "2003MNRAS.338...85M;2005ApJS..161..304M;2009ApJ...696..370W;2009MNRAS.395..342B;2019MNRAS.483..392D",
            Ncl,
        )
    )  # bibcodes for references the data came from
    regions.append(np.repeat("", Ncl))  # galactic subregion designation (if any)
    Vcs.append(np.interp(Rgc, Rgc2, Vrot))  # circular velocity at Rgc
    omegas.append(np.interp(Rgc, Rgc2, omega))  # orbital angular frequency = Vc/Rgc
    kappas.append(np.interp(Rgc, Rgc2, kappa))  # epicyclic frequency
    rho_tidals.append(np.interp(Rgc, Rgc2, rho_tidal))
    DGCs.append(np.repeat(dist_dict[g], Ncl))
    sigma_stars.append(sigma_star)
    sigma_gass.append(sigma_gas)
    sigma_SFRs.append(np.repeat(sigma_SFR_dict[g], Ncl))
    if g == "LMC":
        sigma_1Ds.append(np.repeat(8.0, Ncl))  # Indu 2015A&A...573A.136I
    else:
        sigma_1Ds.append(sigma_1D)  # Indu 2015A&A...573A.136I

# %% [markdown]
# # M83 from Silva-Villa catalogue (Ryon 2015)

# %%
Rgc2, Vrot = np.loadtxt("M83_RGC_vs_Vc.csv").T  # Lundgren 2004
omega = Vrot / Rgc2
kappa = np.sqrt(2 * omega / Rgc2 * np.gradient(Rgc2**2 * omega) / np.gradient(Rgc2))
rho_enc = (
    3 * (4 * omega**2 - kappa**2) / (4 * np.pi) * 2.325e-4
)  # 3 * Menc / (4*np.pi*(1e3*Rgc2)**3)

data = Table.read("Ryon2015_M83.vot", format="votable")
e_Reff = data["e_logReff"]
cut = e_Reff < 1  # no meaningful constraint on Reff
data = data[cut]
gamma = data["Eta"] * 2
# data = data[gamma>2.6] # eliminate shallow profiles with unreliable radii as in Ryon 2015
# gamma = data["Eta"]*2
e_gamma = data["e_Eta"] * 2
e_Reff = data["e_logReff"]
age = data["logAge"]
# data = data[cut]
age = 10 ** data["logAge"]
mass = 10 ** data["logMass"]
mmin = 10 ** data["b_logMass"]
mmax = 10 ** data["B_logMass"]
reff = 10 ** data["logReff"]
error_r = data["e_logReff"]
error_rho = np.sqrt((3 * error_r) ** 2 + (np.log10(mmax / mmin) / 2) ** 2)
Rgc = data["Dgc"]

Rgc3, sigma_gas = np.loadtxt("M83_RGC_vs_Sigmagas.csv").T  # Lundgren 2004
sigma_gas = np.interp(Rgc, Rgc3, sigma_gas)
Rgc3, sigma_SFR = np.loadtxt("M83_RGC_vs_SigmaSFR.csv").T  # Barnes 2014
sigma_SFR = np.interp(Rgc, Rgc3, sigma_SFR)
Rgc3, sigma_star = np.loadtxt("M83_RGC_vs_SigmaStar.csv").T  # Barnes 2014
sigma_star = np.interp(Rgc, Rgc3, sigma_star)
Rgc3, sigma_1D = np.loadtxt("M83_RGC_vs_Sigma1D.csv").T  # Lundgreen 2004
sigma_1D = np.interp(Rgc, Rgc3, sigma_1D)


galaxy_name.append(np.repeat("M83", len(mass)))
references.append(
    np.repeat("2004A&A...422..865L;2014ApJ...789..126B;2015MNRAS.452..525R", len(mass))
)  # bibcodes for references the data came from
masses.append(mass)  # best mass estimate
reffs.append(reff)  # projected effective (half light) radius
reff_error.append(e_Reff)
Rgcs.append(Rgc)  # galactocentric radius
regions.append(np.repeat("", len(mass)))  # galactic subregion designation (if any)
Vcs.append(np.interp(Rgc, Rgc2, Vrot))  # circular velocity at Rgc
omegas.append(np.interp(Rgc, Rgc2, omega))  # orbital angular frequency = Vc/Rgc
kappas.append(np.interp(Rgc, Rgc2, kappa))  # epicyclic frequency
rho_tidals.append(np.interp(Rgc, Rgc2, rho_enc))
ages.append(age)
sigma_gass.append(sigma_gas)
sigma_SFRs.append(sigma_SFR)
sigma_stars.append(sigma_star)
sigma_1Ds.append(sigma_1D)
DGCs.append(np.repeat(4500, len(mass)))  # Madore 2003

# %% [markdown]
# # galaxies from LEGUS catalogue and Brown & Gnedin 2021 fits

# %%
catalog = Table.read("cluster_sizes_brown_gnedin_21.txt", format="ascii.ecsv")

from astropy import units as u
from astropy.coordinates import SkyCoord

galaxy_coords = {
    "ngc5194": (202.469583, 47.2),
    # "ngc628": (24.17393801526, 15.78364097564),
    "ngc4449": (187.04632499999997, 44.093558333333334),
}

galaxy_angles = {
    "ngc5194": 12.0,  # Hu 2013
    # "ngc628": 20.9, # Salak 2019
    "ngc4449": 64,  # Hunter 1999
}
galaxy_incs = {
    "ngc5194": 20.3,  # Hu 2013
    # "ngc628": 8.5, # Walter 2008
    "ngc4449": 45,  # Hunter 1999
}

galaxy_rotationcurve = {
    "ngc5194": "Meidt2013_M51_R_vs_Vrot.csv",
    # "ngc628": "NGC628_RGC_vs_Vrot.csv",
}

# galaxy_vdisp = {
# "ngc5194": "M51_RGC_vs_Sigma1D.csv", # from Pety 2013ApJ...779...43P
# "ngc628": "NGC628_RGC_vs_Sigma1D.csv" # from Mogosi 2016AJ....151...15M
# }

# data = Table.read("Leroy2008.vot",format='votable')
# rad = data["r"][data["Name"]=="NGC 5194"]
# sigma_star = data["Sigma_"][data["Name"]=="NGC 5194"]
# sigma_gas = (data["SigmaH2"]+data["SigmaHI"])[data["Name"]=="NGC 5194"]
# sigma_SFR = (data["FUV_24"])[data["Name"]=="NGC 5194"]
# np.savetxt("M51_RGC_vs_SigmaGas.csv",np.c_[rad,sigma_gas],header="# From Leroy  2008AJ....136.2782L")
# np.savetxt("M51_RGC_vs_SigmaStar.csv",np.c_[rad,sigma_star],header="# From Leroy  2008AJ....136.2782L")
# np.savetxt("M51_RGC_vs_SigmaSFR.csv",np.c_[rad,sigma_SFR/1e4],header="# From Leroy  2008AJ....136.2782L")

# rad = data["r"][data["Name"]=="NGC  628"]
# sigma_star = data["Sigma_"][data["Name"]=="NGC 628"]
# sigma_gas = (data["SigmaH2"]+data["SigmaHI"])[data["Name"]=="NGC  628"]
# sigma_SFR = (data["FUV_24"])[data["Name"]=="NGC 628"]
# np.savetxt("NGC628_RGC_vs_SigmaGas.csv",np.c_[rad,sigma_gas],header="# From Leroy  2008AJ....136.2782L")
# np.savetxt("NGC628_RGC_vs_SigmaStar.csv",np.c_[rad,sigma_star],header="# From Leroy  2008AJ....136.2782L")
# np.savetxt("NGC628_RGC_vs_SigmaSFR.csv",np.c_[rad,sigma_SFR/1e4],header="# From Leroy  2008AJ....136.2782L")

# galaxy_sigma_gas = {"ngc5194": "M51_RGC_vs_SigmaGas.csv",} #"ngc628": "NGC628_RGC_vs_SigmaGas.csv"}
# galaxy_sigma_SFR = {"ngc5194": "M51_RGC_vs_SigmaSFR.csv",} #"ngc628": "NGC628_RGC_vs_SigmaSFR.csv"}
# galaxy_sigma_star = {"ngc5194": "M51_RGC_vs_SigmaStar.csv"}#, "ngc628": "NGC628_RGC_vs_SigmaStar.csv"}

PHANGS_galaxies = ["NGC628", "NGC1433", "NGC1566", "NGC3351", "NGC7793"]

for g in np.unique(catalog["galaxy"]):
    name = g.upper()  # .replace("NGC")
    if not include_galaxies_without_env_data:
        if not (g in galaxy_coords.keys() or name in PHANGS_galaxies):
            continue
    # print(name, name in PHANGS_galaxies)
    if name == "NGC5194":
        name = "M51"
    # parse the LEGUS mass errors
    catalog["mass_msun_e-"] = catalog["mass_msun"] - catalog["mass_msun_min"]
    catalog["mass_msun_e+"] = catalog["mass_msun_max"] - catalog["mass_msun"]

    mask = (
        catalog["reliable_radius"] & catalog["reliable_mass"] & (catalog["galaxy"] == g)
    )  # & (catalog["power_law_slope"] > 1.3) #& (catalog["mass_msun"] > 5e3)

    subset = catalog[mask]
    Ncl = len(subset)
    print(name, Ncl)

    if g in galaxy_incs.keys():
        # get the clusters with reliable radii and masses.
        X0, Y0 = galaxy_coords[g]
        X, Y = catalog["RA"], catalog["Dec"]
        Rgc = np.array(
            [
                correct_rgc(
                    ICRS("%gd" % i, "%gd" % j),
                    glx_ctr=ICRS("%gd" % X0, "%gd" % Y0),
                    glx_PA=Angle(galaxy_angles[g], unit=u.deg),
                    glx_incl=Angle(galaxy_incs[g], unit=u.deg),
                    glx_dist=Distance(
                        catalog["galaxy_distance_mpc"][0] * 1e3, unit=u.kpc
                    ),
                ).value
                for i, j in zip(X, Y)
            ]
        )

        Rgc = Rgc[mask]
    elif name in PHANGS_galaxies:
        # print("Getting PHANGS coordinates for" + name)
        data = Table.read("PHANGS/" + name + "_hexagon_1p5kpc.ecsv")
        X0, Y0 = data.meta["RA_DEG"], data.meta["DEC_DEG"]
        X, Y = catalog["RA"], catalog["Dec"]
        PA = data.meta["PA_DEG"]
        INC = data.meta["INCL_DEG"]
        Rgc = np.array(
            [
                correct_rgc(
                    ICRS("%gd" % i, "%gd" % j),
                    glx_ctr=ICRS("%gd" % X0, "%gd" % Y0),
                    glx_PA=Angle(PA, unit=u.deg),
                    glx_incl=Angle(INC, unit=u.deg),
                    glx_dist=Distance(
                        catalog["galaxy_distance_mpc"][0] * 1e3, unit=u.kpc
                    ),
                ).value
                for i, j in zip(X, Y)
            ]
        )
        Rgc = Rgc[mask]
    else:
        Rgc = np.repeat(np.nan, Ncl)

    try:
        if name in PHANGS_galaxies:
            data = Table.read("PHANGS/" + name + "_annulus_0p5kpc.ecsv")
            Rgc2 = data["r_gal"]
            Vrot = data["V_circ_CO21_URC"]
            if name == "NGC7793":
                Rgc2, Vrot = np.loadtxt(name + "_RGC_vs_Vrot.csv").T
        else:
            Rgc2, Vrot = np.loadtxt(name + "_RGC_vs_Vrot.csv").T
        omega = Vrot / Rgc2
        kappa = np.sqrt(
            2 * omega / Rgc2 * np.gradient(Rgc2**2 * omega, Rgc2)
        )  # / np.gradient(Rgc2))
        rho_enc = 3 * (4 * omega**2 - kappa**2) / (4 * np.pi) * 2.325e-4
        Vrot = np.interp(Rgc, Rgc2, Vrot)
        omega = np.interp(Rgc, Rgc2, omega)
        kappa = np.interp(Rgc, Rgc2, kappa)
        rho_enc = np.interp(Rgc, Rgc2, rho_enc)
    except:
        print("Rotation curve not found for " + name)
        Vrot = omega = kappa = rho_enc = np.repeat(np.nan, Ncl)
        # Rgc2, Vrot = np.array([0,100]), np.array([np.nan,np.nan])

    try:
        if name in PHANGS_galaxies:
            # print("opening " + "PHANGS/"+name+"_annulus_0p5kpc.ecsv")
            data = Table.read("PHANGS/" + name + "_annulus_0p5kpc.ecsv")

            Rgc3 = data["r_gal"]
            sigma_mol = data["Sigma_mol"]
            sigma_atom = data["Sigma_atom"]
            sigma_mol[np.isnan(sigma_mol)] = 0.0
            sigma_atom[np.isnan(sigma_atom)] = 0.0
            sigma_gas = sigma_mol + sigma_atom
            sigma_gas = np.interp(Rgc, Rgc3, sigma_gas)
        else:
            Rgc3, sigma_gas = np.loadtxt(name + "_RGC_vs_SigmaGas.csv").T
            sigma_gas = np.interp(Rgc, Rgc3, sigma_gas)
    except:
        sigma_gas = np.repeat(np.nan, Ncl)

    try:
        if name in PHANGS_galaxies:
            data = Table.read("PHANGS/" + name + "_annulus_0p5kpc.ecsv")
            Rgc3 = data["r_gal"]
            sigma_SFR = data["Sigma_SFR_HaW4recal"]
            if np.all(np.isnan(sigma_SFR)):
                sigma_SFR = data["Sigma_SFR_FUVW4"]
        else:
            Rgc3, sigma_SFR = np.loadtxt(name + "_RGC_vs_SigmaSFR.csv").T
        sigma_SFR = np.interp(Rgc, Rgc3, sigma_SFR)
    except:
        sigma_SFR = np.repeat(np.nan, Ncl)

    try:
        if name in PHANGS_galaxies:
            data = Table.read("PHANGS/" + name + "_annulus_0p5kpc.ecsv")
            Rgc3 = data["r_gal"]
            sigma_star = data["Sigma_star"]
        else:
            Rgc3, sigma_star = np.loadtxt(name + "_RGC_vs_SigmaStar.csv").T
        sigma_star = np.interp(Rgc, Rgc3, sigma_star)
    except:
        sigma_star = np.repeat(np.nan, Ncl)

    try:
        if name in PHANGS_galaxies:
            data = Table.read("PHANGS/" + name + "_annulus_0p5kpc.ecsv")
            Rgc3 = data["r_gal"]
            sigma_1D = data["<vdisp_mol_pix_150pc>"]
        else:
            Rgc3, sigma_1D = np.loadtxt(
                name + "_RGC_vs_Sigma1D.csv"
            ).T  # np.linspace(0.1,Rgc.max(),100)
    except:
        print(
            "Velocity dispersion not found for "
            + name
            + "; assuming 11km/s per Leroy 2008"
        )
        Rgc3, sigma_1D = np.array([0, 100]), np.array([11.0, 11.0])
    sigma_1D = np.interp(Rgc, Rgc3, sigma_1D)

    mcl, rcl = subset["mass_msun"], subset["r_eff_pc"]  # /0.75
    e_logm = np.log10(subset["mass_msun_max"] / subset["mass_msun_min"]) / 2
    e_logr = np.log10(subset["r_eff_pc_e+"] / subset["r_eff_pc_e-"]) / 2
    error_rho = np.sqrt((3 * e_logr) ** 2 + e_logm**2)
    error_rho = np.sqrt(error_rho**2)  # + 0.3**2)
    rhocl = 3 * mcl / (8 * np.pi * rcl**3)

    mass = mcl
    reff = rcl
    age = subset["age_yr"]
    galaxy_name.append(np.repeat(name, len(mass)))
    references.append(
        np.repeat(
            "2008AJ....136.2782L;2023ApJ...945L..19S;2021MNRAS.508.5935B", len(mass)
        )
    )  # bibcodes for references the data came from
    masses.append(mass)  # best mass estimate
    reffs.append(reff)  # projected effective (half light) radius
    reff_error.append(e_logr)
    Rgcs.append(Rgc)  # galactocentric radius
    regions.append(np.repeat("", len(mass)))  # galactic subregion designation (if any)
    Vcs.append(Vrot)  # circular velocity at Rgc
    omegas.append(omega)  # orbital angular frequency = Vc/Rgc
    kappas.append(kappa)  # epicyclic frequency
    rho_tidals.append(rho_enc)
    ages.append(age)
    sigma_gass.append(sigma_gas)
    sigma_SFRs.append(sigma_SFR)
    sigma_stars.append(sigma_star)
    sigma_1Ds.append(sigma_1D)
    DGCs.append(
        subset["galaxy_distance_mpc"] * 1e3
    )  # TRGB distances from Sabbi 2018, except NGC1566

# %%
# for x in galaxy_name, references, masses, reffs, Rgcs, regions, Vcs, omegas, kappas, rho_tidals:#, sigma_gass, sigma_SFRs:
# x = np.concatenate(x)
print(np.unique(np.concatenate(galaxy_name)))
compilation = Table()
compilation["GalaxyName"] = np.concatenate(galaxy_name)
compilation["GalaxyName"] = np.array([str.upper(g) for g in compilation["GalaxyName"]])
compilation["References"] = np.concatenate(references)
compilation["Mass"] = np.concatenate(masses)
compilation["Reff"] = np.concatenate(reffs)
compilation["error_logReff"] = np.concatenate(reff_error)
compilation["Rgc"] = np.concatenate(Rgcs)
compilation["Region"] = np.concatenate(regions)
compilation["Vcirc"] = np.concatenate(Vcs)
compilation["Omega"] = np.concatenate(omegas)
compilation["Age"] = np.concatenate(ages)
compilation["Kappa"] = np.concatenate(kappas)
compilation["RhoTidal"] = np.concatenate(rho_tidals)
compilation["SigmaGas"] = np.concatenate(sigma_gass)
compilation["SigmaStar"] = np.concatenate(sigma_stars)
compilation["SigmaSFR"] = np.concatenate(sigma_SFRs)
compilation["GalaxyDistance"] = np.concatenate(DGCs)
compilation["SigmaGas1D"] = np.concatenate(sigma_1Ds)

from astropy.io import ascii

ascii.write(compilation, "Cluster_Compilation.dat", format="basic", overwrite=True)
