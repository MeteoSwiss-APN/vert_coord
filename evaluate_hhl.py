from turtle import color
from matplotlib import markers
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import psyplot.project as psy
import iconarray
import click
from ipdb import set_trace
import sys
from pathlib import Path
import os

# example commands:
# python evaluate_hhl.py --print_info
# python evaluate_hhl.py --print_info --model cosmo-1 --file /store/s83/swester/grids/const_modinterim.nc


def get_min_max(arr):
    """Minimum and maximum.

    Args:
        arr (np array): array

    """
    minimum = np.min(arr)
    maximum = np.max(arr)

    return minimum, maximum


def ind_from_latlon(lats, lons, lat, lon, verbose=False):
    """Find the nearest neighbouring index to given location.

    Args:
        lats (1d array):            Latitude grid
        lons (1d array):            Longitude grid
        lat (float):                Latitude of location
        lon (float):                Longitude of location
        verbose (bool, optional):   Print information. Defaults to False.

    Returns:
        int     Index of nearest grid point.
    """
    dist = [
        np.sqrt((lats[i] - lat) ** 2 + (lons[i] - lon) ** 2) for i in range(len(lats))
    ]
    ind = np.where(dist == np.min(dist))[0][0]

    if verbose:
        print(f"Closest ind: {ind}")
        print(f" Given lat: {lat:.3f} vs found lat: {lats[ind]:.3f}")
        print(f" Given lat: {lon:.3f} vs found lon: {lons[ind]:.3f}")

    return ind


def get_poi(lats, lons):
    """Points of interest for analysis.

    Args:
        lats
        lons

    Returns:
        pd dataframe

    """
    poi = pd.DataFrame(
        columns=["mtblanc", "zrh", "pay", "visp", "ulr", "sav", "duf", "cic"],
        index=["long_name", "ind", "h_real"],
    )

    # indeces of specific locations
    ind_mtblanc = ind_from_latlon(lats, lons, 45.83267, 6.86437)
    ind_zrh = ind_from_latlon(lats, lons, 47.46218, 8.54458)
    ind_pay = ind_from_latlon(lats, lons, 46.81291, 6.94418)
    ind_visp = ind_from_latlon(lats, lons, 46.29861, 7.88004)
    ind_ulr = ind_from_latlon(lats, lons, 46.50568, 8.30610)
    ind_sav = ind_from_latlon(lats, lons, 44.276917, 8.546750)
    ind_duf = ind_from_latlon(lats, lons, 45.93692, 7.86675)
    ind_cic = ind_from_latlon(lats, lons, 45.72350, 8.61444)

    poi["mtblanc"].long_name = "Mt Blanc"
    poi["zrh"].long_name = "Zürich"
    poi["pay"].long_name = "Payerne"
    poi["visp"].long_name = "Visp"
    poi["ulr"].long_name = "Ulrichen"
    poi["sav"].long_name = "Savona"
    poi["duf"].long_name = "Dufourspitze"
    poi["cic"].long_name = "Cicognola"

    poi["mtblanc"].ind = ind_mtblanc
    poi["zrh"].ind = ind_zrh
    poi["pay"].ind = ind_pay
    poi["visp"].ind = ind_visp
    poi["ulr"].ind = ind_ulr
    poi["sav"].ind = ind_sav
    poi["duf"].ind = ind_duf
    poi["cic"].ind = ind_cic

    poi["mtblanc"].h_real = 4808.0
    poi["zrh"].h_real = 422.0
    poi["pay"].h_real = 491.0
    poi["visp"].h_real = 646.0
    poi["ulr"].h_real = 1345.0
    poi["sav"].h_real = 0.0
    poi["duf"].h_real = 4634.0
    poi["cic"].h_real = 197.0

    return poi


def n_sum_up_to(dz, top):
    """Number of dz until sum of dz reaches top.

    Args:
        dz (1d array): level thicknesses
        top (upper limit): upper limit

    Return:
        int: Number of dz

    """
    cumsum = np.cumsum(dz)
    return np.sum(cumsum < top)


def info_minmax(hsurf):
    """Print information elevation of heighest and lowest grid cell.

    Args:
        hsurf (2d array): surface height
    """

    # HSURF
    min_hsurf, max_hsurf = get_min_max(hsurf)
    print(f"Lowest point in domain: {min_hsurf:.2f} m asl.")
    print(f"Highest point in domain: {max_hsurf:.2f} m asl.")

    return


def info_hhl(hhl, poi):
    """Print information about height of half levels to screen.

    Args:
        hhl (3d array): height of half levels (including surface)
        poi (pd dataframe): points of interest
    """

    # print info on highest and lowest elevation grid cell
    info_minmax(hhl[-1, :])

    # Number of levels
    nlev = hhl.shape[0]

    for location in poi:
        p = poi[location]
        print(f"\n{p.long_name}:")
        print(f"   Real elevation: {int(p.h_real)} m asl.")
        print(f"   Model elevation: {int(hhl[-1, p.ind])} m asl.")
        for i in range(1, 12):
            print(f"   {i}. level: {hhl[-i,p.ind]:.2f} m asl")
        print(f"   ...")
        for i in range(2, -1, -1):
            print(f"   {nlev-i}. level: {hhl[i,p.ind]:.2f} m asl")

    return


def info_dz(hsurf, poi, dz):
    """Print information about level distribution to screen.

    Args:
        hsurf (2d array): surface height
        poi (pd dataframe): point of interest
        dz (3d array): level thickness
    """

    # print info on highest and lowest elevation grid cell
    info_minmax(hsurf)

    for location in poi:
        p = poi[location]
        print(f"\n{p.long_name}:")
        print(f"   Real elevation: {int(p.h_real)} m asl.")
        print(f"   Model elevation: {int(hsurf[p.ind])} m asl.")
        for i in range(1, 11):
            print(f"   {i}. dz: {dz[-i,p.ind]:.2f} m")
        nl_500 = n_sum_up_to(dz[::-1, p.ind], 500)
        print(f"   Levels in lowest 500 m above ground: {nl_500}")

    # lowest level
    dz_low = dz[-1, :]
    print(f"\nLowest level:")
    print(f"  - maximum thickness: {max(dz_low):.2f}m")
    print(f"  - minimum thickness: {min(dz_low):.2f}m")
    set_trace()
    # max_dz_neighbours(hsurf)

    print(f"")


def mapplot_coord_surf(file, grid_file, lev):

    print(f"Plotting surface elevation of {file}.")
    print(f"Corresponding grid file: {grid_file}.")

    ds = iconarray.combine_grid_information(file, grid_file)

    # get nlev, find name
    try:
        nlev = ds.HHL.values.shape[0]
        name_height = "HHL"
    except AttributeError:
        nlev = ds.HEIGHT.values.shape[0]
        name_height = "HEIGHT"
    zz = nlev - lev

    surf_map = ds.psy.plot.mapplot(
        name=name_height,
        xgrid=None,
        ygrid=None,
        # bounds=np.arange(0, 101, 5),
        # norm=norm,
        cticksize="small",
        # cmap=cmap,
        borders=True,
        lakes=True,
        rivers=True,
        projection="robin",
        # map_extent=[5.5, 11.0, 45.5, 48.0],
        title=f"Elevation of {lev}. coordinate surface",
        clabel="m asl",
        z=zz,  # specify specific level
    )

    user = os.getlogin()
    out_dir = Path(f"/scratch/{user}/vert_coord/figures")
    out_dir.mkdir(exist_ok=True)
    out_name = Path(out_dir, f"altitude_{lev}_coordinate_surface.png")
    plt.savefig(out_name)
    print(f"Saved as: {out_name}")


def plot_dz(dz, poi, model, exp, out_dir):
    """Plot dz vs altitude.

    Args:
        dz (np.array): level thickness
        poi (pd.Dataframe): points of interest
        model (str): model name
        exp (str): identifier for specific configuration
        out_dir (str): output directory

    """
    # plot "nominal" dz, i.e. at location of savona
    ii = poi["sav"].ind
    dz_nominal = dz[::-1, ii]
    n_lev = len(dz_nominal)

    # cut into lower and upper region
    cut = 22
    dz_lower = dz_nominal[:cut]
    dz_upper = dz_nominal[(cut - 1) :]

    ### dz vs level-index ###
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 6))

    # plot all levels
    ax1.plot(np.arange(1, cut + 1), dz_lower, marker="o", color="purple")
    ax1.plot(np.arange(cut, n_lev + 1), dz_upper, marker="o", color="yellowgreen")
    ax1.set_xlabel("# Level")
    ax1.set_ylabel("dz [m]")
    ax1.set_title("All levels")
    ax1.set_ylim(0, 1200)
    ax1.set_xlim(0, n_lev)

    # plot only lower atmosphere
    ax2.plot(np.arange(1, cut + 1), dz_lower, marker="o", color="purple")
    ax2.set_xlabel("# Level")
    ax2.set_title("Lower atmosphere")
    ax2.set_ylim(min(dz_lower) - 2, max(dz_lower) + 5)
    ax2.set_xlim(0, cut)

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_name = Path(out_dir, f"dz_vs_level_{model}_{exp}.png")
    plt.tight_layout()
    plt.savefig(out_name)
    print(f"Saved as: {out_name}")

    plt.clf()

    ### altitude vs dz ###
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 8))

    # calculate cumulative sum of dz
    dz_cumsum = np.cumsum(dz_nominal)

    # plot all levels
    ax1.plot(dz_lower, dz_cumsum[:cut], marker="o", color="purple")
    ax1.plot(dz_upper, dz_cumsum[(cut - 1) :], marker="o", color="yellowgreen")
    ax1.set_xlabel("dz [m]")
    ax1.set_ylabel("Altitude [m asl]")
    ax1.set_title("All levels")
    ax1.set_xlim(0, 1200)
    ax1.set_ylim(0, 20500)

    # plot only lower atmosphere
    ax2.plot(dz_lower, dz_cumsum[:cut], marker="o", color="purple")
    ax2.set_xlabel("dz [m]")
    ax2.set_title("Lower atmosphere")
    ax2.set_xlim(min(dz_lower) - 2, max(dz_lower) + 2)
    ax2.set_ylim(0, dz_cumsum[(cut - 1)])

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_name = Path(out_dir, f"altitude_vs_dz_{model}_{exp}.png")
    plt.tight_layout()
    plt.savefig(out_name)
    print(f"Saved as: {out_name}")


@click.command()
@click.option(
    "--file",
    default="/store/mch/msopr/swester/grids/HEIGHT_ICON-1E.nc",
    help="Netcdf file containing HHL to be analysed.",
    type=str,
)
@click.option(
    "--grid_file",
    default="/store/s83/tsm/ICON_INPUT/icon-1e_dev/ICON-1E_DOM01.nc",
    help="Netcdf file containing grid information.",
    type=str,
)
@click.option(
    "--model",
    default="icon-1",
    help="originating model",
    type=str,
)
@click.option(
    "--config",
    default="ref",
    help="specific configuration",
    type=str,
)
@click.option(
    "--print_dz",
    is_flag=True,
    default=False,
)
@click.option(
    "--print_hhl",
    is_flag=True,
    default=False,
)
@click.option(
    "--plot_surf",
    is_flag=True,
    default=False,
)
@click.option(
    "--lev",
    type=int,
    default=1,
)
@click.option(
    "--loc",
    type=str,
    multiple=True,
    default=["all"],
)
@click.option(
    "--create_plots",
    is_flag=True,
    default=False,
)
@click.option(
    "--out_dir",
    default="/scratch/swester/vert_coord/figures",
    type=str,
)
def evaluate_hhl(
    file,
    grid_file,
    model,
    config,
    print_dz,
    print_hhl,
    plot_surf,
    lev,
    loc,
    create_plots,
    out_dir,
):

    print(f"Evaluating {file}")

    # load file, retrieve relevant variables
    ds = psy.open_dataset(file).squeeze()

    if "icon" in model:
        try:
            hhl = ds.HEIGHT.values
            lats = ds.clat_1.values
            lons = ds.clon_1.values
        except AttributeError:
            hhl = ds.HHL.values
            lats = ds.clat.values
            lons = ds.clon.values

    if max(lats) < 2:
        lats = np.rad2deg(lats)
        lons = np.rad2deg(lons)

    elif "cosmo" in model:
        hhl_3d = ds.HEIGHT.values
        s1 = hhl_3d.shape[1]
        s2 = hhl_3d.shape[2]
        hhl = hhl_3d.reshape(hhl_3d.shape[0], s1 * s2)
        lats = ds.lat_1.values.flatten()
        lons = ds.lon_1.values.flatten()

    else:
        print("Model type unknown!")
        sys.exit(1)

    hsurf = hhl[-1, :]
    dz = hhl[:-1, :] - hhl[1:, :]
    hfl = hhl[1:, :] + dz[:, :] / 2

    # load points of interest
    all_poi = get_poi(lats, lons)
    if loc[0] == "all":
        poi = all_poi
    else:
        poi = all_poi[list(loc)]

    if print_dz:
        info_dz(hsurf, poi, dz)

    if print_hhl:
        info_hhl(hhl, poi)

    if plot_surf:
        if "icon" in model:
            mapplot_coord_surf(file, grid_file, poi, lev)
        else:
            print(f"No mapplot available for {model}.")

    if create_plots:
        plot_dz(dz, poi, model, config, out_dir)


if __name__ == "__main__":
    evaluate_hhl()
