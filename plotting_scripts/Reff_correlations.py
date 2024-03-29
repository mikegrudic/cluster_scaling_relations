from sys import argv
import numpy as np
import matplotlib.pyplot as plt

# from fix_log_ticks import fix_log_ticks
from astropy.table import Table
from scipy.special import erfinv

p_to_sigma = lambda p: (erfinv(1 - p) * 2**0.5).clip(0, 8)

from multiprocessing import Pool
from scipy.stats import spearmanr


data = Table.read("../Cluster_Compilation.dat", format="ascii.basic")

PHANGS_LEGUS = ["NGC628", "NGC1433", "NGC1566", "NGC3351", "NGC7793"]

ages = data["Age"]
max_age = 3e8
min_age = 1e7
cut = np.bool_((ages < max_age) * (ages > min_age) * (data["Reff"] > 0))
data = data[cut]
galaxies = data["GalaxyName"]


completeness_limits = {
    "M31": 2e3,
    "M51": 5e3,
    "LMC": 1e2,
    "SMC": 1e2,
    "M83": 1e4,
    "NGC4449": 5e3,
    "NGC628": 5e3,
    "M33": 2e3,
    "M82": 1.5e4,
}
legus_galaxies = [
    "NGC45",
    "NGC4395",
    "NGC1433",
    "NGC1705",
    "M51",
    "NGC1313",
    "NGC1566",
    "NGC3344",
    "NGC3738",
    "NGC4449",
    "NGC4656",
    "NGC5474",
    "NGC628",
    "NGC6503",
    "NGC7793",
    "ugc1249",
    "ugc7408",
    "NGC3351",
    "NGC4242",
]
for g in galaxies:
    if not g in completeness_limits.keys():
        completeness_limits[g] = 5e3

ages = data["Age"]
sigma_gas = data["SigmaGas"]
fGMC = (1 + 0.025 * (sigma_gas / 100) ** -2) ** -1
phi_P_hat = 10 - 8 * fGMC
phi_P = 3
alpha_vir = 1.3
sigma_GMC = np.sqrt(10.0 / 3 * phi_P * phi_P_hat * alpha_vir) * sigma_gas
sigma_GMC2 = 8 + 1.8 * (sigma_gas) ** 0.9
sigma_SFR = data["SigmaSFR"]
sigma_star = data["SigmaStar"]
rho_tidal = (
    3 * (4 * data["Omega"] ** 2 - data["Kappa"] ** 2) / (4 * np.pi) * 2.325e-4
)  # * 125
tidal_field = 4 * data["Omega"] ** 2 - data["Kappa"] ** 2
sigma_1D = data["SigmaGas1D"]
Q = sigma_1D * data["Kappa"] / data["SigmaGas"] * 0.074
m = data["Mass"]
r = data["Reff"]
sigma_eff = m / (2 * np.pi * r**2)
Rgc = data["Rgc"]
rho_eff = 3 * m / (8 * np.pi * r**3 * 1.3**3)
gal = data["GalaxyName"]


gal = data["GalaxyName"]

sigmaz = data["SigmaGas1D"]
# exponential disk scale lengths
scale_lengths = {
    "M31": 5.91,  # Seigar 2008
    "M51": 4.2,  # Beckman 1996 M51 B band value
    "M33": 1.8,  # (Verley et al. 2007; Ferguson et al. 2007; Verleyet al. 2009; Corbelli et al. 2014)
    "M82": 0.8,  # Mayya 2009
    "LMC": 1.6,  # Alves 2000
    "SMC": 0.8,  # Nidever 2011
    "M83": 1.7,  # Barnes 2014
    "NGC628": 3.2,  # Salo 2015
    "NGC4449": 3.3 / 3.2,  # Sacchi 2016
}

scale_lengths = np.array(
    [scale_lengths[g] if g in scale_lengths else np.nan for g in gal]
)
rho_star = data["SigmaStar"] / (0.54 * 1000 * scale_lengths)  # in msun pc^-2
G = 4.3e-3  # msun-pc-km/s units
A = np.sqrt(G * sigma_gas**2 / (2 * rho_star * sigmaz**2))
H = sigmaz / np.sqrt(4 * np.pi * G * rho_star) / (A + np.sqrt(1 + A**2))
rho0 = sigma_gas * (2 * np.pi) ** -0.5 / H
P_midplane = sigmaz**2 * rho0
rho_compressed = P_midplane / (0.2**2 * (sigma_SFR / 1e-2) ** 0.25)
size_linewidth = sigmaz**2 / (2 * H)
k1, k2 = 1.0, 3.0 / 5
sigma_GMC3 = 3 * k1 / (2 * np.pi * k2) * size_linewidth / G
rho_clump = (
    9
    * np.pi
    / 20
    * 2
    * 0.5**2
    * G
    * sigma_GMC3**2
    / (0.2**2 * (data["SigmaSFR"] / 1e-2) ** 0.25)
)
gamma_GMC = sigmaz / (sigma_GMC3 * rho0)

dependent_variables = (r,)  # sigma_eff, rho_eff
independent_variables = (
    m,
    ages,
    Rgc,
    sigma_star,
    sigma_gas,
    sigma_SFR,
    tidal_field,
    sigma_1D,
)  # sigma_GMC3, rho_clump, gamma_GMC, sigma_gas * Q**2
itags = [
    r"$M$",
    r"$\tau$",
    r"$R_{\rm GC}$",
    r"$\Sigma_\star$",
    r"$\Sigma_{\rm g}$",
    r"$\Sigma_{\rm SFR}$",
    r"$4\Omega^2-\kappa^2$",
    r"$\sigma_{\rm z}$",
]  # ,r"$\Sigma_{\rm cloud}$", r"$\rho_{\rm clump}$", r"$\gamma_{\rm GMC}$", r"$\Sigma_{\rm CK21}$"]

fig, ax = plt.subplots(1, 1, figsize=(3, 6))

corrs = []
ps = []
gals = []
Ns = []

do_PHANGS_LEGUS = False
PHANGS_LEGUS = ["NGC628", "NGC1433", "NGC1566", "NGC3351", "NGC7793"]
for g in np.unique(galaxies):
    # if(galaxies==g).sum() < 10: continue
    idx = (galaxies == g) * (m > completeness_limits[g])  # *np.isfinite(Rgc)
    print(g, idx.sum())
    if idx.sum() < 10:
        continue
    gals.append(g)
    Ns.append(idx.sum())
    # if not g=="M82": continue
    for n, i in enumerate(independent_variables):
        x = i[idx]
        if np.all(np.isnan(x)):
            corr = np.nan
            p = np.nan
        else:
            y = r[idx]
            idx2 = np.isfinite(x) * np.isfinite(y)
            x, y = x[idx2], y[idx2]
            corr, p = spearmanr(y, x)
        corrs.append(corr)
        ps.append(p)

if do_PHANGS_LEGUS:
    idx = np.array([g in PHANGS_LEGUS for g in galaxies]) * (
        m > completeness_limits["NGC628"]
    )
    gals.append("LEGUS-PHANGS")
    Ns.append(idx.sum())
    for n, i in enumerate(independent_variables):
        x = i[idx]
        if np.all(np.isnan(x)):
            corr = np.nan
            p = np.nan
        else:
            y = r[idx]
            idx2 = np.isfinite(x) * np.isfinite(y)
            x, y = x[idx2], y[idx2]
            corr, p = spearmanr(y, x)
        corrs.append(corr)
        ps.append(p)


corrs = np.array(corrs).reshape(len(gals), len(independent_variables))
ps = np.array(ps).reshape(len(gals), len(independent_variables))
order = np.array(Ns).argsort()[::-1]
corrs = corrs[order]
ps = ps[order]
gals = np.array(gals)[order]
Ns = np.array(Ns)[order]

corrs = np.ma.array(corrs, mask=np.isnan(corrs))
# corrs[ps > 0.05] = np.nan
cmap = plt.cm.RdBu
cmap.set_bad("#DDDDDD", 1.0)
# corrs[ps>0.05] = 1.
plot = ax.matshow(corrs, cmap=cmap, vmin=-1, vmax=1)
for i in range(corrs.shape[0]):
    for j in range(corrs.shape[1]):
        print(i, j, ps[i, j])
        if np.isfinite(ps[i, j]):
            ax.text(
                j,
                i,
                (r"$%2.2g\sigma$" % p_to_sigma(ps[i, j])).replace("0.", "."),
                fontsize=5,
                horizontalalignment="center",
                verticalalignment="center",
                color=("black" if ps[i, j] < 0.05 else "grey"),
            )

ax.minorticks_off()
ax.tick_params(
    axis="x", which="both", bottom=False, top=False, labelbottom=False, labeltop=True
)
ax.tick_params(axis="y", which="both", left=False, right=False)
from mpl_toolkits.axes_grid1 import make_axes_locatable

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("bottom", "5%", pad="0%")
cbar = fig.colorbar(plot, pad=0, cax=cax, orientation="horizontal")
cbar.ax.tick_params(labelsize=8)
cbar.set_label(r"Correlation with $R_{\rm eff}$", size=8)
ax.set_yticks(range(len(gals)))
ax.grid(which="minor")
ax.set_yticklabels(
    [
        (
            r"$\mathbf{%s}$" % str.upper(g) + r" ($N=%d$)" % n
            if (g in PHANGS_LEGUS or g == "LEGUS-PHANGS")
            else str.upper(g) + r" ($N=%d$)" % n
        )
        for g, n in zip(gals, Ns)
    ],
    fontsize=6,
)
ax.set_xticklabels([""] + itags, rotation=90, fontsize=8)
plt.savefig("Reff_correlations.pdf", bbox_inches="tight")
