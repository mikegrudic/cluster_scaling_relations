from matplotlib import pyplot as plt
import numpy as np
from astropy.table import Table
from pandas import read_csv
import palettable

datasets = [
    # "MvdM05 (Local Group)",
    # "Ryon+15 (M83)",
    # "Ryon+17 (NGC 628,1313)",
    "BG21 (LEGUS)",
    # "Gatto+21 (SMC)",
]
DATA_PATH = "../data/"


def get_gamma_values(ref, agemin=1e7, agemax=3e8):
    """Returns the EFF slope from a given dataset"""
    agecut = lambda age: (age > agemin) * (age < agemax)

    if ref == "Ryon+15 (M83)":
        data = Table.read(DATA_PATH + "Ryon2015_M83.vot", format="votable")
        gamma = data["Eta"] * 2
        e_gamma = data["e_Eta"] * 2
        age = 10 ** data["logAge"]
        cut = agecut(age)
        return [["M83", gamma[cut], e_gamma[cut]]]
    if ref == "Ryon+17 (NGC 628,1313)":
        result = []
        for gal in "NGC628", "NGC1313":
            data = Table.read(DATA_PATH + f"Ryon2017_{gal}.vot", format="votable")
            gamma = data["eta"] * 2
            e_gamma = data["e_eta"] * 2
            age = 10 ** data["logAge"]
            cut = agecut(age)
            result.append([gal, gamma[cut], e_gamma[cut]])
        return result
    if ref == "BG21 (LEGUS)":
        data = Table.read(
            DATA_PATH + "cluster_sizes_brown_gnedin_21.txt", format="ascii.ecsv"
        )
        result = []
        for gal in np.unique(data["galaxy"]):
            gamma = 2 * data["power_law_slope"]
            e_gamma = data["power_law_slope_e-"] + data["power_law_slope_e+"]
            age = data["age_yr"]
            cut = agecut(age) * (data["galaxy"] == gal)
            result.append([gal, gamma[cut], e_gamma[cut]])
        return result
    if ref == "Gatto+21 (SMC)":
        gal = "SMC"
        data = read_csv(DATA_PATH + "gatto_2021_SMC_table1.csv")
        data2 = read_csv(DATA_PATH + "gatto_2021_SMC_tableC1.csv")
        gamma = data2["gamma"]
        e_gamma = data2["gamma_std"]
        age = 10 ** data["logt"]
        cut = agecut(age)
        return [["SMC", gamma[cut], e_gamma[cut]]]
    if ref == "MvdM05 (Local Group)":
        result = []
        t1 = Table.read("../data/MM05_cluster_params.vot")
        t2 = Table.read("../data/MM05_cluster_photometry.vot")
        # t4 = Table.read("fits?J%2FApJS%2F161%2F304%2Ftable14.dat")

        galaxies = [n.split("-")[0] for n in t1["Cluster"]]

        for g in np.unique(galaxies):
            model = t1["Mod"]
            mask = (model == "PL") * (
                np.array(galaxies) == g
            )  # * np.isfinite(t1["Rh"])  # select EFF fits

            gamma = t1["W0_gamma"][mask] - 1
            e_gamma = t1["e_W0_gamma"][mask]
            agedict = dict(zip(t2["Cluster"], 10 ** t2["logAge"]))
            age = np.array([agedict[n] for n in t1["Cluster"][mask]])
            cut = agecut(age)
            if cut.sum():
                result.append([g, gamma, e_gamma])
        # print(result)
        return result


fig, ax = plt.subplots(figsize=(4, 4))
ax.set_prop_cycle("color", palettable.cartocolors.qualitative.Bold_5.mpl_colors)
gamma_bins = np.linspace(1, 5, 21)
for d in datasets:
    data = get_gamma_values(d)
    gammas = []
    for gal, gamma, e_gamma in data:
        gamma = np.array(gamma)
        # if len(gamma) < 30:
        # continue
        gamma = gamma[np.isfinite(gamma)]
        gammas.append(gamma)
    gammas = np.concatenate(gammas)
    ax.ecdf(gammas, label=d)  # gal.upper() + f" ({d})")

ax.legend(labelspacing=0)
ax.set(xlim=[1, 5], xlabel=r"$\gamma$", ylabel=r"Fraction $<\gamma$")
plt.savefig("Gamma_EDF.png", dpi=800)
plt.clf()


fig, ax = plt.subplots(figsize=(4, 4))
ax.set_prop_cycle("color", palettable.cartocolors.qualitative.Bold_5.mpl_colors)
gamma_bins = np.linspace(1, 6, 16)
for d in datasets:
    data = get_gamma_values(d)
    gammas = []
    for gal, gamma, e_gamma in data:
        gamma = np.array(gamma)
        gamma = gamma[np.isfinite(gamma)]
        gammas.append(gamma)
    gammas = np.concatenate(gammas)
    ax.hist(
        gammas,
        gamma_bins,
        label=d,
        density=True,
        color="black",
        histtype="step",
    )  # gal.upper() + f" ({d})")

ax.legend(labelspacing=0)
ax.set(xlim=[1, 6], xlabel=r"$\gamma$", ylabel=r"Fraction $<\gamma$")
plt.savefig("Gamma_PDF.png", dpi=800)
