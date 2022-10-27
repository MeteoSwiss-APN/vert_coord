# construct vertical coordinate from surface upwards
# Author: Stephanie Westerhuis
# Date: October 2022
####################################################

# python packages
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import psyplot.project as psy
import iconarray
import click
from ipdb import set_trace
import sys
from pathlib import Path
import os

# home-made functions

# example commands:
#
# python construct_hhl.py --file /store/s83/swester/vert_coord_files/icon-1-alps/const_galchen.nc


def calc_vct_ab(n_levels, top_height, stretch_fac, h_flat):
    """Calculate vct_a and vct_b.

    FORTRAN CODE:
    DO jk = 1, nlevp1

        z_flat = flat_height  ! from 'sleve_nml'

        x1 = REAL( nlev + 1 - jk, wp) / REAL( nlev, wp)
        z_height  = top_height * x1 * ( stretch_fac * x1 + 1.0_wp-stretch_fac )
        vct_a(jk) = z_height
        IF ( z_height >= z_flat) THEN
        vct_b(jk) = 0.0_wp
        ELSE
        vct_b(jk) = (z_flat - z_height)/z_flat
        ENDIF

    ENDDO

    Args:
        n_levels (int):         Number of levels, including surface, usually 81.
        top_height (float):     Altitude of model domain top.
        stretch_fac (float):    Stretching factor.
        h_flat (float):         Altitude above which levels are flat.

    Returns:
        2 1-dimensional vectors: vct_a and vct_b
    """

    vec_n_levels = np.arange(1, n_levels + 1)
    x1 = (vec_n_levels - 1) / (n_levels - 1)
    vct_a = top_height * x1 * (stretch_fac * x1 + 1 - stretch_fac)
    vct_b = (h_flat - vct_a) / h_flat

    # find lowest level which is above h_flat
    i_flat = np.min(vec_n_levels[vct_a > h_flat]) - 1

    # fill vct_b from this level onwards with 0s
    vct_b[i_flat:] = 0

    # obtain the same index order as in the FORTRAN code
    vct_a = np.flip(vct_a)
    vct_b = np.flip(vct_b)
    return vct_a, vct_b


@click.command()
@click.option(
    "--file",
    default="/store/s83/swester/vert_coord_files/icon-1-alps/external_parameter_icon_alps_R19B08_mch.nc",
    help="REQUIRED: Netcdf file containing orography.",
    type=str,
)
@click.option(
    "--grid_file",
    default="/store/s83/swester/vert_coord_files/icon-1-alps/alps_DOM01.nc",
    help="REQUIRED: Netcdf file containing grid information.",
    type=str,
)
@click.option(
    "--n_levels",
    default=81,
    help="Number of vertical full levels (including surface)",
    type=int,
)
@click.option(
    "--h_flat",
    default=16000.0,
    help="Altitude [m] above which levels are flat.",
    type=float,
)
@click.option(
    "--top_height",
    default=22000.0,
    help="Model domain top altitude [m].",
    type=float,
)
@click.option(
    "--stretch_fac",
    default=0.65,
    help="Stretch factor, between 0 and 1.",
    type=int,
)
@click.option(
    "--type_vct_a",
    default="2nd-order",
    help="Type of nominal level distribution.",
    type=str,
)
@click.option(
    "--type_vct_b",
    default="linear",
    help="Type of decay of orographic signal",
    type=str,
)
@click.option(
    "--out_dir",
    help="Change output directory for figures.",
    default="/scratch/swester/vert_coord/hhl/",
    type=str,
)
@click.option(
    "--verify",
    help="Print information for verification to screen.",
    is_flag=True,
    default=False,
    type=bool,
)
def construct_hhl(
    file,
    grid_file,
    n_levels,
    h_flat,
    top_height,
    stretch_fac,
    type_vct_a,
    type_vct_b,
    out_dir,
    verify,
):

    print(f"\n--- Constructing HHL for {file}\n")
    print(f"--- Specified settings:")
    print(f"    n_levels:       {n_levels}")
    print(f"    h_flat:         {h_flat} m")
    print(f"    top_height:     {top_height} m")
    print(f"    stretch_fac:    {stretch_fac}")
    print(f"    type_vct_a:     {type_vct_a}")
    print(f"    type_vct_b:     {type_vct_b}")

    # load file, retrieve relevant variables
    ds = xr.open_dataset(file).squeeze()

    try:
        hsurf = ds.HSURF.values
        lats = ds.clat_1.values
        lons = ds.clon_1.values
    except AttributeError:
        try:
            hsurf = ds.HSURF.values
            lats = ds.clat.values
            lons = ds.clon.values
        except AttributeError:
            try:
                hsurf = ds.topography_c.values
                lats = ds.clat.values
                lons = ds.clon.values
            except AttributeError:
                print("--- names for height field, latitudes or longitudes unknown!")
                print(f"--- check: {file}")
                sys.exit(0)

    if max(lats) < 2:
        lats = np.rad2deg(lats)
        lons = np.rad2deg(lons)

    n_cells = len(hsurf)
    print(f"\n--- This domain comprises {n_cells} cells.")

    # calculate vct_a and vct_b
    vct_a, vct_b = calc_vct_ab(n_levels, top_height, stretch_fac, h_flat)
    print("\nvct_a:")
    print(vct_a)
    print("\nvct_b:")
    print(vct_b)

    print("\n Papperlapapp.")


if __name__ == "__main__":
    construct_hhl()
