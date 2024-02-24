from matplotlib import pyplot as plt
from astropy.table import Table
import numpy as np
from scipy.stats import theilslopes

fig, ax = plt.subplots(figsize=(4, 4))

datasets = {
    "EFF": "Cluster_Compilation_EFFcorr.dat",
    "Truncated EFF": "Cluster_Compilation.dat",
}
for label, path in datasets.items():
    data = Table.read(path, format="ascii.basic")
    cut = (
        data["References"]
        == "2008AJ....136.2782L;2023ApJ...945L..19S;2021MNRAS.508.5935B"
    )
    gal = data["GalaxyName"]
    # cut *= gal == "M51"
    data = data[cut]
    M_solar = data["Mass"]
    Reff_pc = data["Reff"]

    ax.loglog(
        M_solar,
        Reff_pc,
        ".",
        marker={"EFF": "o", "Truncated EFF": "s"}[label],
        label=label,
        markersize=1,
        alpha=0.7,
    )

    cut = (M_solar > 1e4) * (Reff_pc > 0.1) * (Reff_pc < 100)
    print(theilslopes(np.log10(M_solar[cut]), np.log10(Reff_pc[cut])))
ax.legend()
ax.set(xlim=[1e3, 1e7], ylim=[0.1, 100])
plt.show()
