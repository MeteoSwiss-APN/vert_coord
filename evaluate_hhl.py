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
from utils import open_ds_regular
from utils import retrieve_vars_print
from utils import retrieve_vars_icon_regular
from utils import retrieve_vars_cosmo_regular
from printing import info_minmax
from printing import info_hhl
from printing import info_dz
from printing import info_max_dzdc
from plotting import transect_hhl
from plotting import transect_topo
from plotting import transect_topo_regular
from plotting import mapplot_coord_surf
from plotting import profile_dz

# COSMO-1:
# python evaluate_hhl.py --print_dz --model cosmo-1
#   --file /store/s83/swester/grids/const_modinterim.nc


@click.command()
@click.option(
    "--model",
    default="icon",
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
    help="OPTIONAL: Second file for comparison.",
    type=str,
)
@click.option(
    "--grid_file",
    default="/scratch/e1000/meteoswiss/scratch/swester/input_icon/grids/icon-1e_dev/ICON-1E_DOM01.nc",
    help="REQUIRED FOR ICON: Netcdf file containing grid information.",
    type=str,
)
@click.option(
    "--config",
    default="ref",
    help="SPECIFY: Name of configuration.",
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
    help="OPTIONAL: Select specific location(s) only, defaults to all. Select from 'mtblanc', 'zrh', 'pay', 'visp', 'ulr', 'sav', 'duf', 'cic'",
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
    help="SPECIFY: Index of level from buttom upwards, surface = 1.",
    type=int,
    default=1,
)
@click.option(
    "--radius",
    help="SPECIFY: Number of grid cells around location to be included in plot.",
    type=int,
    default=30,
)
@click.option(
    "--vmin",
    help="SPECIFY: Minimum elevation.",
    type=float,
    default=0.0,
)
@click.option(
    "--vmax",
    help="SPECIFY: Maximum elevation.",
    type=float,
    default=4500.0,
)
@click.option(
    "--plot_hhl",
    help="SELECT: Plot vertical transect of vertical coordinate surfaces.",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "--plot_topo",
    help="SELECT: Plot topography.",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "--plot_ddz",
    help="SELECT: Plot d delta_z / dz.",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "--out_dir",
    help="SPECIFY: Change output directory for figures.",
    default="figures",
    type=str,
)
@click.option(
    "--verify",
    help="SELECT: Print information for verification to screen.",
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
    plot_ddz,
    out_dir,
    verify,
):

    print(f"\nEvaluating {file}\n")

    out_dir = parse_out_dir(out_dir)

    if print_dz:
        lats, lons, hhl, hsurf, dz = retrieve_vars_print(file, model)
        poi = get_poi(loc, lats, lons)
        print("Printing dz...\n")
        info_dz(hsurf, poi, dz, lev)

    # if print_min_dz:
    #    lats, lons, hhl, hsurf, dz = retrieve_vars_print(file, model)
    #    poi = get_poi(loc, lats, lons)
    #    print("Printing min dz...\n")
    #    info_min_dz(hsurf, poi, dz, lev)

    if print_hhl:
        lats, lons, hhl, hsurf, dz = retrieve_vars_print(file, model)
        poi = get_poi(loc, lats, lons)
        print("Printing hhl...\n")
        info_hhl(hhl, poi, lev)

    if print_max_dzdc:

        if "cosmo" in model:
            print("Not implemented for COSMO.")
            sys.exit()

        lats, lons, hhl, hsurf, dz = retrieve_vars_print(file, model)
        poi = get_poi(loc, lats, lons)

        print("Printing maximum elevation difference...\n")
        info_max_dzdc(hhl, grid_file, poi, lev, lats, lons, verify)

    if plot_surf:
        print("Plotting vertical coordinate surface...\n")
        if "icon" in model:
            ds, ds_grid = open_icon_ds(file, grid_file)
            mapplot_coord_surf(ds, config, out_dir, lev, loc, radius, vmin, vmax)
        else:
            print(f"No mapplot available for {model}.")

    if plot_hhl:
        if "icon" in model:
            ds, ds_grid = open_icon_ds(file, grid_file)
            transect_hhl(ds, ds_grid, loc, config, out_dir, lev)
        else:
            print(f"No mapplot available for {model}.")

    #    transect_hhl(hhl, neighbour_ind, poi, config, out_dir, lev)

    if plot_ddz:
        lats, lons, hhl, hsurf, dz = retrieve_vars_print(file, model)
        poi = get_poi(loc, lats, lons)
        print("Plotting d delta_z / dz...\n")
        profile_dz(dz, hhl, poi, loc, config, out_dir)

    if plot_topo:
        if "regular" in model:
            ds = open_ds_regular(file)
            lats, lons, hsurf = retrieve_vars_icon_regular(ds)
            poi = get_poi(loc, lats, lons, model="icon_regular")
            transect_topo_regular(lats, lons, hsurf, poi, radius, config, out_dir)
        elif "icon" in model:
            print("first re-integrate")
            sys.exit()
            neighbour_ind = ds_grid.neighbor_cell_index.values
            transect_topo(hhl, neighbour_ind, ds, poi, config, out_dir, lev)
        elif "cosmo" in model:
            ds = open_ds_regular(file)
            lats, lons, hsurf = retrieve_vars_cosmo_regular(ds)
            poi = get_poi(loc, lats, lons, model=model)
            transect_topo_regular(lats, lons, hsurf, poi, radius, config, out_dir)
        else:
            print("Not sure which model to take.")


if __name__ == "__main__":
    evaluate_hhl()
