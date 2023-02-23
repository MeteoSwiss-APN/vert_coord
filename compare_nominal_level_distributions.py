# python packages
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from ipdb import set_trace
from pathlib import Path

# homemade packages
from utils import get_poi
from utils import retrieve_vars_print
from printing import info_hhl, info_dz

# some generic settings
out_dir = "/scratch/e1000/meteoswiss/scratch/swester/vert_coord/figures/nominal_levels"
c1 = "skyblue"
c2 = "greenyellow"
c3 = "darkblue"
c4 = "forestgreen"
locations = ["sav", "mtblanc"]
max_levels = 60
max_levels_zoom = 20

# files
f1 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/dwd20.nc"
f2 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/cosmo20.nc"
f3 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/dwd15.nc"
f4 = "/store/mch/msopr/swester/vert_coord_files/compare_nominal_levels/cosmo15.nc"

# get hhl and points of interest
print("=== retrieving lats, lons, hhl, dz")
lats, lons, hhl1, hsurf1, dz1 = retrieve_vars_print(f1, model="icon")
lats, lons, hhl2, hsurf2, dz2 = retrieve_vars_print(f2, model="icon")
lats, lons, hhl3, hsurf3, dz3 = retrieve_vars_print(f3, model="icon")
lats, lons, hhl4, hsurf4, dz4 = retrieve_vars_print(f4, model="icon")

print("=== loading points of interest")
poi = get_poi(loc=locations, lats=lats, lons=lons)

p = poi["mtblanc"]
for lev in range(1, max_levels_zoom):
    print(f"=== Thickness of {lev}. level:")
    print(
        f"    DWD20: {dz1[-lev,p.ind]:.2f}m, COSMO20: {dz2[-lev,p.ind]:.2f}m, DWD15: {dz3[-lev,p.ind]:.2f}m, COSMO15: {dz4[-lev,p.ind]:.2f}m"
    )


levels = np.arange(81, 0, -1)
# plotting
for location in locations:

    p = poi[location]

    # altitude vs N level
    fig, ax = plt.subplots(figsize=(4, 6))
    ax.plot(
        levels[-max_levels:],
        hhl1[-max_levels:, p.ind],
        color=c1,
        marker="o",
        markersize=2,
        label="107: DWD_20m",
    )
    ax.plot(
        levels[-max_levels:],
        hhl2[-max_levels:, p.ind],
        color=c2,
        marker="o",
        markersize=2,
        label="108: COSMO_20m",
    )
    ax.plot(
        levels[-max_levels:],
        hhl3[-max_levels:, p.ind],
        color=c3,
        marker="o",
        markersize=2,
        label="109: DWD_15m",
    )
    ax.plot(
        levels[-max_levels:],
        hhl4[-max_levels:, p.ind],
        color=c4,
        marker="o",
        markersize=2,
        label="110: COSMO_15m",
    )
    ax.set_xlabel("Level")
    ax.set_ylabel("Altitude [masl]")
    ax.set_xlim(0, lev)
    ax.set_ylim(hhl1[-1, p.ind], hhl1[-max_levels, p.ind])
    ax.legend()
    ax.set_title("Altitude of model levels")
    out_name = Path(out_dir, f"altitude_levels_{location}")
    plt.savefig(out_name, dpi=200, bbox_inches="tight")
    print(f"Saved as: {out_name}.png")

    ax.set_xlim(0, max_levels_zoom)
    ax.set_ylim(hhl1[-1, p.ind], hhl1[-max_levels_zoom, p.ind])
    out_name = Path(out_dir, f"altitude_levels_{location}_zoom")
    plt.savefig(out_name, dpi=200, bbox_inches="tight")
    print(f"Saved as: {out_name}.png")

    plt.close()

    # N level vs dz
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(
        dz1[-max_levels:, p.ind],
        levels[-max_levels:],
        color=c1,
        marker="o",
        markersize=2,
        label="107: DWD_20m",
    )
    ax.plot(
        dz2[-max_levels:, p.ind],
        levels[-max_levels:],
        color=c2,
        marker="o",
        markersize=2,
        label="108: COSMO_20m",
    )
    ax.plot(
        dz3[-max_levels:, p.ind],
        levels[-max_levels:],
        color=c3,
        marker="o",
        markersize=2,
        label="109: DWD_15m",
    )
    ax.plot(
        dz2[-max_levels:, p.ind],
        levels[-max_levels:],
        color=c4,
        marker="o",
        markersize=2,
        label="110: COSMO_15m",
    )
    ax.set_xlabel("Level thickness [m]")
    ax.set_ylabel("N level")
    ax.set_xlim(0, dz1[-max_levels, p.ind])
    ax.set_ylim(0, max_levels)
    ax.legend()
    ax.set_title("Thickness of model levels")
    out_name = Path(out_dir, f"dz_level_{location}")
    plt.savefig(out_name, dpi=200, bbox_inches="tight")
    print(f"Saved as: {out_name}.png")

    ax.set_xlim(0, dz1[-max_levels_zoom, p.ind])
    ax.set_ylim(0, max_levels_zoom)
    out_name = Path(out_dir, f"dz_level_{location}_zoom")
    plt.savefig(out_name, dpi=200, bbox_inches="tight")
    print(f"Saved as: {out_name}.png")

    plt.close()
