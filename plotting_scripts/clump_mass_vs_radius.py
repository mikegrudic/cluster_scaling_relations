import numpy as np
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Dark2_5

# from fix_log_ticks import fix_log_ticks

fig, ax = plt.subplots(1, 1, sharey=True, figsize=(4, 4))
ax.set_prop_cycle("color", Dark2_5.mpl_colors)
# ax2 = ax.twinx()
###### OBSERVATIONAL POINTS ##############################################################################
from astropy import table
import matplotlib.pyplot as plt
import numpy as np

catalog = table.Table.read("../data/cluster_sizes_brown_gnedin_21.txt", format="ascii.ecsv")

# # Conditions for the Brick: M=1.3e5, R=2.9 (Longmore 2012), T = 27K
# sigma_brick = (1.3e5 * 0.5 / (np.pi * (2.9)**2))
# rhocl_brick = 200 *(sigma_brick/50)**2 / (100/10)
# #ax.text(5e2,rhocl_brick*0.7,r"  $\,\,\,\,\,\,\,\rho_{\rm clump}$" + "\n" +  "The Brick",fontsize=8) #r"$T=50\rm K$; $\Sigma_{\rm GMC}=3000M_\odot\rm pc^{-2}$",fontsize=5)
# ax.text(5e2,rhocl_brick*0.6,r"$T=100\rm K;$" + "\n" + r"$\Sigma_{\rm GMC}=3000M_\odot\rm pc^{-2}$",fontsize=8)
# plot_ellipse(1e3,rhocl_brick,1,1)

# # ATLASGAL clumps ###############################################################

from astropy.io.votable import parse_single_table

table = parse_single_table("../clumpdata/Urquhart_2018_atlasgal_clumps.vot")
# data = table
stage = table.array["EvolType"]
# for s in np.unique(stage):
# cut =
Mcl = 10 ** table.array["logMclump"]
Rcl = (
    table.array["Rad"] * 1.538
)  # factor 1.538 corrects for conversion between half-mass and effective radius ("sigma" of the Gaussian)
RGC = table.array["RGC"]
rhocl = Mcl * 3 / (8 * np.pi * Rcl**3)
ax.scatter(Mcl, Rcl, s=0.1, label="Urquhart 2018 (dust)")

print(np.polyfit(np.log10(Mcl), np.log10(Rcl), 1))

# LMC 12CO clumps

GMC_names = "30Dor", "GMC1", "GMC104", "N59C", "A439", "PCC"

mleaves = []
rleaves = []
mbranch = []
rbranch = []
for n in "12", "13":
    # mleaves = []
    # rleaves = []
    # mbranch = []
    # rbranch = []
    for name in GMC_names:
        leaves = np.int_(np.loadtxt("../clumpdata/" + name + "_%s_leaves.txt" % n))
        branches = np.int_(np.loadtxt("../clumpdata/" + name + "_%s_branches.txt" % n))
        print("../clumpdata/" + name + "_%s_physprop_add.txt" % n)
        data = np.loadtxt("../clumpdata/" + name + "_%s_physprop_add.txt" % n)
        mleaves.append(data[leaves, 8])
        rleaves.append(data[leaves, 2])
        mbranch.append(data[branches, 8])
        rbranch.append(data[branches, 2])
mleaves = np.concatenate(mleaves)
rleaves = np.concatenate(rleaves)
mbranch = np.concatenate(mbranch)
rbranch = np.concatenate(rbranch)
# ax.scatter(mleaves,rleaves,s=1,label=r"Wong 2019 LMC $^{%s}$CO (Leaf)"%n,marker='d')
# ax.scatter(mbranch,rbranch,s=1,label=r"Wong 2019 LMC $^{%s}$CO (Branch)"%n,marker='*')
ax.scatter(
    np.concatenate([mleaves, mbranch]),
    np.concatenate([rleaves, rbranch]),
    s=1,
    label=r"Wong 2019 LMC (CO)",
    marker="d",
)
# ax.scatter(mbranch,rbranch,s=1,label=r"Wong 2019 LMC $^{%s}$CO (Branch)"%n,marker='*')

# FKM points
m, r = np.loadtxt("../clumpdata/FKM2010_clumps_red.csv").T
ax.scatter(m, r, s=2, marker="s", label="Shirley 2003 (CS)")
m, r = np.loadtxt("../clumpdata/FKM2010_clumps_black.csv").T
ax.scatter(m, r, s=2, marker="o", label=r"Fontani 2005 (C$^{17}$O, dust)")
m, rho = np.loadtxt("../clumpdata/Faundez_2004_clumps.txt").T
r = (3 * m / (4 * np.pi * rho)) ** (1.0 / 3)
r *= 10**0.4 / r.max()  # correcting to definition of radius used in Faundez
ax.scatter(m, r, marker="^", label="Faundez 2004 (dust)", s=4)


for nh in 1e2, 1e4, 1e6:
    ms = np.array([10, 3e5])
    rho = nh / 26
    r = (3 * ms / (8 * np.pi * rho)) ** (1.0 / 3)
    ax.plot(ms, r, color="black", ls="dashed", lw=1, zorder=100)
    ax.text(
        ms[1] * 0.1,
        r[1] * 0.1**0.33 * 0.87,
        r"$\mathbf{log \, n=%d}$" % int(np.log10(nh) + 0.5),
        fontsize=8,
        rotation=32,
        ha="left",
    )

for sigma in 1e2, 1e4:
    ms = np.array([10, 3e5])
    r = (ms / (np.pi * sigma)) ** (1.0 / 2)
    ax.plot(ms, r, color="black", ls="dashdot", lw=1, zorder=100)
    if sigma == 1e4:
        ax.text(
            ms[1] * 0.1,
            r[1] * 0.1**0.5 * 1.2,
            r"$\mathbf{log\, \Sigma=%d}$" % int(np.log10(sigma) + 0.5),
            fontsize=8,
            rotation=43,
            ha="left",
        )
    else:
        ax.text(
            ms[1] * 0.07,
            r[1] * 0.07**0.5 * 1.2,
            r"$\mathbf{log\,\Sigma=%d}$" % int(np.log10(sigma) + 0.5),
            fontsize=8,
            rotation=43,
            ha="left",
        )

ms = np.array([1e3, 3e5])
ax.fill_between(
    ms,
    2.5 * (ms / 1e4) ** 0.28 * 10**-0.24,
    2.5 * (ms / 1e4) ** 0.28 * 10**0.24,
    color="black",
    label="LEGUS Clusters $\pm \sigma$",
    alpha=0.3,
    zorder=-1000,
)
# ax.fill_between(ms,(ms/1e4)**0.28*10**-0.24,color='grey')


ax.set(
    xscale="log",
    yscale="log",
    xlim=[10, 3e5],
    ylim=[0.1, 20],
    xlabel=r"Mass ($M_\odot$)",
    ylabel="Radius (pc)",
)
ax.legend(labelspacing=0, fontsize=8)  # ,ncol=2)
plt.savefig("clump_mass_vs_radius.png", bbox_inches="tight", dpi=400)
