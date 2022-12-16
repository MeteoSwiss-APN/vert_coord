# evaluate vertical coordinate from ICON simulation
# Author: Stephanie Westerhuis
# Date: September 2022
###################################################

# python packages
import numpy as np
import pandas as pd
import xarray as xr
import click
from ipdb import set_trace
import sys
from pathlib import Path
import os

# home-made functions
from utils import get_min_max
from utils import ind_from_latlon
from utils import get_poi
from utils import n_sum_up_to
from utils import parse_out_dir
from utils import open_icon_ds
from plotting import transect_hhl
from plotting import transect_topo
from plotting import mapplot_coord_surf

# COSMO-1:
# python evaluate_hhl.py --print_dz --model cosmo-1
#   --file /store/s83/swester/grids/const_modinterim.nc


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


def info_hhl(hhl, poi, lev):
    """Print information about height of half levels to screen.

    Args:
        hhl (3d array): height of half levels (including surface)
        poi (pd dataframe): points of interest
        lev (int): N of level counted from surface upwards
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
        for i in range(1, lev):
            print(f"   {i}. level: {hhl[-i,p.ind]:.2f} m asl")
        print(f"   ...")
        for i in range(2, -1, -1):
            print(f"   {nlev-i}. level: {hhl[i,p.ind]:.2f} m asl")

    return


def info_dz(hsurf, poi, dz, lev):
    """Print information about level distribution to screen.

    Args:
        hsurf (2d array): surface height
        poi (pd dataframe): point of interest
        dz (3d array): level thickness
        lev (int): N of level (from surface upwards)
    """

    # print info on highest and lowest elevation grid cell
    info_minmax(hsurf)

    for location in poi:
        p = poi[location]
        print(f"\n{p.long_name}:")
        print(f"   Real elevation: {int(p.h_real)} m asl.")
        print(f"   Model elevation: {int(hsurf[p.ind])} m asl.")
        for i in range(1, lev):
            print(f"   {i}. dz: {dz[-i,p.ind]:.2f} m")
        nl_500 = n_sum_up_to(dz[::-1, p.ind], 500)
        print(f"   Levels in lowest 500 m above ground: {nl_500}")

    # lowest level
    dz_low = dz[-1, :]
    print(f"\nLowest level:")
    print(f"  - maximum thickness: {max(dz_low):.2f}m")
    print(f"  - minimum thickness: {min(dz_low):.2f}m")
    print(f"")


def info_max_dzdc(hhl, grid_file, poi, lev, lats, lons, verify=False):
    """Print information about maximum z-difference between cells.

    Args:
        hhl (2d np.array):      height of half levels
        grid_file (str):        path to grid file
        poi (pd.DataFrame):     points of interest
        lev (int):              N of level (upwards)
        lats (1d np.array):     latitude
        lons (1d np.array):     longitude
    """

    # load grid file
    grid = iconarray.open_dataset(grid_file)

    # neighbour indeces (3 per cell) in fortran start
    #  at 1, hence subtract 1 to make pythonic
    neighs = grid.neighbor_cell_index.values - 1

    # specified coordinate surface
    surf = hhl[-lev, :]

    if verify:
        ii = 3789
        print(f"--- Value of cell {ii} on level {lev}:")
        print(f"    {surf[ii]:.2f}")

    # for simpler code reading:
    # split neighbour vectors for 1st/2nd/3rd neighbour
    # n0 is just a vector with index as value
    n0 = np.arange(0, len(surf))
    n1 = neighs[0, :]
    n2 = neighs[1, :]
    n3 = neighs[2, :]

    if verify:
        print(f"--- Neighbouring cell indeces:")
        print(f"    {n1[ii]}, {n2[ii]}, {n3[ii]}")
        print(f"--- Neighbouring cell values:")
        print(f"    {surf[n1[ii]]:.2f}, {surf[n2[ii]]:.2f}, {surf[n3[ii]]:.2f}")

    # fill "no-neighbour" cells with itself index (n0) such that
    # calculation of difference to neighbouring cell
    # will access itself and dz = 0

    # a) retrieve indices
    to_fill_n1 = np.where(n1 < 0)[0]
    to_fill_n2 = np.where(n2 < 0)[0]
    to_fill_n3 = np.where(n3 < 0)[0]

    if verify:
        print(f"--- Will fill {len(to_fill_n1)} cells indicating 1st neighbour.")
        print(f"--- Will fill {len(to_fill_n2)} cells indicating 2nd neighbour.")
        print(f"--- Will fill {len(to_fill_n3)} cells indicating 3rd neighbour.")

    # b) fill them
    n1[to_fill_n1] = n0[to_fill_n1]
    n2[to_fill_n2] = n0[to_fill_n2]
    n3[to_fill_n3] = n0[to_fill_n3]

    # create 3 new fields with same shape as "surf"
    # and put value of 1st/2nd/3rd neighbour of cell into cell
    surf_n1 = surf[n1]
    surf_n2 = surf[n2]
    surf_n3 = surf[n3]

    if verify:
        print(f"--- Surface of neighbour 1 at cell {ii}:")
        print(f"    {surf_n1[ii]:.2f}")
        print(f"--- Surface of neighbour 2 at cell {ii}:")
        print(f"    {surf_n2[ii]:.2f}")
        print(f"--- Surface of neighbour 3 at cell {ii}:")
        print(f"    {surf_n3[ii]:.2f}")

    # calculate absolute difference fields
    dz_n1 = np.abs(surf - surf_n1)
    dz_n2 = np.abs(surf - surf_n2)
    dz_n3 = np.abs(surf - surf_n3)

    if verify:
        print(f"--- dz to neighbour 1 at cell {ii}:")
        print(f"    {dz_n1[ii]}")
        print(f"--- dz to neighbour 2 at cell {ii}:")
        print(f"    {dz_n2[ii]}")
        print(f"--- dz to neighbour 3 at cell {ii}:")
        print(f"    {dz_n3[ii]}")

    # finally, determine maximum dz between adjacent cells
    dzdc = np.maximum.reduce([dz_n1, dz_n2, dz_n3])
    max_dzdc = np.max(dzdc)
    max_ii = np.where(dzdc == max_dzdc)[0][0]

    if verify:
        print(f"--- maximum dz from cell {ii} to neighbours:")
        print(f"    {max_dzdc[ii]}")

    print("\n***********************************************************")
    print(f"--- Maximum dz between 2 adjacent grid cells on level {lev}:")
    print(f"    {max_dzdc:.2f}")
    print(f"--- location:")
    print(f"    {lats[max_ii]:.3f}, {lons[max_ii]:.3f}")
    print("\n***********************************************************\n")
    for location in poi:
        p = poi[location]
        print(f"--- Maximum dz to any of the 3 neighbouring cells at {p.long_name}:")
        print(f"    {dzdc[p.ind]:.2f}")
    print("***********************************************************")


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
    "--model",
    default="icon-1",
    help="REQUIRED: Originating model, e.g. 'cosmo' or 'icon'.",
    type=str,
)
@click.option(
    "--file",
    default="/store/mch/msopr/swester/grids/HEIGHT_ICON-1E.nc",
    help="REQUIRED: Netcdf file containing HHL to be analysed.",
    type=str,
)
@click.option(
    "--file2",
    help="Second file for comparison.",
    type=str,
)
@click.option(
    "--grid_file",
    default="/store/s83/swester/vert_coord_files/icon-1-alps/alps_DOM01.nc",
    help="REQUIRED FOR ICON: Netcdf file containing grid information.",
    type=str,
)
@click.option(
    "--config",
    default="ref",
    help="Name to specify configuration.",
    type=str,
)
@click.option(
    "--print_dz",
    help="SELECT: Prints details on levels' thickness.",
    is_flag=True,
    default=False,
)
@click.option(
    "--print_hhl",
    help="SELECT: Prints details on altitude of half levels.",
    is_flag=True,
    default=False,
)
@click.option(
    "--print_max_dzdc",
    help="SELECT: Print maximum elevation difference between adjacent grid cells.",
    is_flag=True,
    default=False,
)
@click.option(
    "--loc",
    help="Print information for specific location(s) only, defaults to all. Select from 'mtblanc', 'zrh', 'pay', 'visp', 'ulr', 'sav', 'duf', 'cic'",
    type=str,
    multiple=True,
    default=["all"],
)
@click.option(
    "--plot_surf",
    help="SELECT: Plots surface elevation in 2d mapplot.",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "--lev",
    help="Specify (number of) level(s), surface = 1.",
    type=int,
    default=1,
)
@click.option(
    "--radius",
    help="Number of grid cells around location to be included in plot.",
    type=int,
    default=30,
)
@click.option(
    "--vmin",
    help="Minimum elevation.",
    type=float,
    default=0.0,
)
@click.option(
    "--vmax",
    help="Maximum elevation.",
    type=float,
    default=4500.0,
)
@click.option(
    "--plot_hhl",
    help="Plot vertical transect of vertical coordinate surfaces.",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "--plot_topo",
    help="Plot topography.",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "--out_dir",
    help="Change output directory for figures.",
    default="vc_scratch",
    type=str,
)
@click.option(
    "--verify",
    help="Print information for verification to screen.",
    is_flag=True,
    default=False,
    type=bool,
)
def evaluate_hhl(
    file,
    file2,
    grid_file,
    model,
    config,
    print_dz,
    print_hhl,
    print_max_dzdc,
    plot_surf,
    lev,
    radius,
    vmin,
    vmax,
    loc,
    plot_hhl,
    plot_topo,
    out_dir,
    verify,
):

    print(f"\nEvaluating {file}\n")

    out_dir = parse_out_dir(out_dir)
    print(out_dir)

    # routines require either the files, xarray datasets,
    #  or numpy arrays of the retrieved variables

    if print_dz:
        lats, lons, hhl, hsurf, dz = retrieve_vars(file, model)
        poi = get_poi(lats, lons, loc)
        print("Printing dz...\n")
        info_dz(hsurf, poi, dz, lev)

    if print_hhl:
        lats, lons, hhl, hsurf, dz = retrieve_vars(file, model)
        poi = get_poi(lats, lons, loc)
        print("Printing hhl...\n")
        info_hhl(hhl, poi, lev)

    if print_max_dzdc:

        if "cosmo" in model:
            print("Not implemented for COSMO.")
            sys.exit()

        lats, lons, hhl, hsurf, dz = retrieve_vars(file, model)
        poi = get_poi(lats, lons, loc)

        print("Printing maximum elevation difference...\n")
        info_max_dzdc(hhl, grid_file, poi, lev, lats, lons, verify)

    if plot_surf:
        print("Plotting vertical coordinate surface...\n")
        if "icon" in model:
            ds = open_icon_ds(file, grid_file)
            mapplot_coord_surf(ds, config, out_dir, lev, loc, radius, vmin, vmax)
        else:
            print(f"No mapplot available for {model}.")

    # if plot_hhl:
    #    transect_hhl(hhl, neighbour_ind, poi, config, out_dir, lev)

    # if plot_topo:
    #  neighbour_ind = ds_grid.neighbor_cell_index.values
    #    transect_topo(hhl, neighbour_ind, ds, poi, config, out_dir, lev)


if __name__ == "__main__":
    evaluate_hhl()
