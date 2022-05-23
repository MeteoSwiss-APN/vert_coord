import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import psyplot.project as psy
import click
from ipdb import set_trace
import sys


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
        columns=["mtblanc", "zrh", "pay", "visp", "ulr", "sav"],
        index=["long_name", "ind", "h_real"],
    )

    # indeces of specific locations
    ind_mtblanc = ind_from_latlon(lats, lons, 45.83267, 6.86437)
    ind_zrh = ind_from_latlon(lats, lons, 47.46218, 8.54458)
    ind_pay = ind_from_latlon(lats, lons, 46.81291, 6.94418)
    ind_visp = ind_from_latlon(lats, lons, 46.29861, 7.88004)
    ind_ulr = ind_from_latlon(lats, lons, 46.50568, 8.30610)
    ind_sav = ind_from_latlon(lats, lons, 44.276917, 8.546750)

    poi["mtblanc"].long_name = "Mt Blanc"
    poi["zrh"].long_name = "Zürich"
    poi["pay"].long_name = "Payerne"
    poi["visp"].long_name = "Visp"
    poi["ulr"].long_name = "Ulrichen"
    poi["sav"].long_name = "Savona"

    poi["mtblanc"].ind = ind_mtblanc
    poi["zrh"].ind = ind_zrh
    poi["pay"].ind = ind_pay
    poi["visp"].ind = ind_visp
    poi["ulr"].ind = ind_ulr
    poi["sav"].ind = ind_sav

    poi["mtblanc"].h_real = 4808.0
    poi["zrh"].h_real = 422.0
    poi["pay"].h_real = 491.0
    poi["visp"].h_real = 646.0
    poi["ulr"].h_real = 1345.0
    poi["sav"].h_real = 0.0

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


def info_dz(hsurf, poi, dz):
    """Print information about level distribution to screen.

    Args:
        hsurf (2d array): surface height
        poi (pd dataframe): point of interest
        dz (3d array): level thickness
    """

    # HSURF
    min_hsurf, max_hsurf = get_min_max(hsurf)
    print(f"Lowest point in domain: {min_hsurf:.2f}")
    print(f"Highest point in domain: {max_hsurf:.2f}")

    for location in poi:
        p = poi[location]
        print(f"\n{p.long_name}:")
        print(f"   Real elevation: {int(p.h_real)} m asl.")
        print(f"   Model elevation: {int(hsurf[p.ind])} m asl.")
        for i in range(1, 11):
            print(f"   {i}. dz: {dz[-i,p.ind]}")
        nl_500 = n_sum_up_to(dz[::-1, p.ind], 500)
        print(f"   Levels in lowest 500m above ground: {nl_500}")

    set_trace()
    # max_dz_neighbours(hsurf)

    print(f"")


@click.command()
@click.option(
    "--file",
    default="/store/mch/msopr/swester/grids/HEIGHT_ICON-1E.nc",
    help="Netcdf file containing HHL to be analysed.",
    type=str,
)
@click.option(
    "--model",
    default="icon",
    help="originating model",
    type=str,
)
@click.option(
    "--print_info",
    is_flag=True,
    default=False,
)
@click.option(
    "--plot_dz",
    is_flag=True,
    default=False,
)
@click.option(
    "--out_dir",
    default="/scratch/swester/vert_coord/",
    type=str,
)
def evaluate_hhl(file, model, print_info, plot_dz, out_dir):

    print(f"Evaluating {file}")

    # load file, retrieve relevant variables
    ds = psy.open_dataset(file).squeeze()

    if model == "icon":
        hhl = ds.HEIGHT.values
        lats = ds.clat_1.values
        lons = ds.clon_1.values

    elif model == "cosmo":
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
    poi = get_poi(lats, lons)

    if print_info:
        info_dz(hsurf, poi, dz)


if __name__ == "__main__":
    evaluate_hhl()
