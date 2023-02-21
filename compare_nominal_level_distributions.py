# python packages
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from ipdb import set_trace
from pathlib import Path

# homemade packages
from utils import get_poi
from utils import retrieve_lats_lons_hhl_icon
from printing import info_hhl

# some generic settings
out_dir = "/scratch/e1000/meteoswiss/scratch/swester/vert_coord/figures/nominal_levels"
c1 = "orchid"
c2 = "slategray"
c3 = "goldenrod"
locations = [
    "sav",
]

# files
f1 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/dwd.nc"
f2 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/cosmo20.nc"
f3 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/cosmo15.nc"

# load netcdf's
print("=== opening datasets")
n1 = xr.open_dataset(f1).squeeze()
n2 = xr.open_dataset(f2).squeeze()
n3 = xr.open_dataset(f3).squeeze()

# get hhl and points of interest
print("=== retrieving lats, lons, hhl")
lats, lons, hhl_1 = retrieve_lats_lons_hhl_icon(n1)
lats, lons, hhl_2 = retrieve_lats_lons_hhl_icon(n2)
lats, lons, hhl_3 = retrieve_lats_lons_hhl_icon(n3)

print("=== loading points of interest")
poi = get_poi(loc=locations, lats=lats, lons=lons)

for hhl in [hhl_1, hhl_2, hhl_3]:
    info_hhl(hhl, poi, lev=15)

levels = np.arange(81, 0, -1)
for location in locations:

    p = poi[location]

    # altitude vs N level
    fig, ax = plt.subplots(figsize=(4, 6))
    ax.plot(
        levels,
        hhl_1[:, p.ind],
        color=c1,
        marker="o",
        markersize=2,
        label="107: DWD",
    )
    ax.plot(
        levels,
        hhl_2[:, p.ind],
        color=c2,
        marker="o",
        markersize=2,
        label="108: COSMO_20m",
    )
    ax.plot(
        levels,
        hhl_3[:, p.ind],
        color=c3,
        marker="o",
        markersize=2,
        label="109: 'COSMO'_15m",
    )
    ax.set_xlabel("Level")
    ax.set_ylabel("Altitude [masl]")
    ax.set_xlim(0, 55)
    ax.set_ylim(0, 6000)
    ax.legend()
    ax.set_title("Altitude of model levels")
    out_name = Path(out_dir, f"altitude_levels_{location}")
    plt.savefig(out_name, dpi=200, bbox_inches="tight")
    print(f"Saved as: {out_name}.png")

    ax.set_xlim(0, 20)
    ax.set_ylim(0, 800)
    out_name = Path(out_dir, f"altitude_levels_{location}_zoom")
    plt.savefig(out_name, dpi=200, bbox_inches="tight")
    print(f"Saved as: {out_name}.png")

    plt.close()
