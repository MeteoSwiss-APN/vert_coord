# utilities to evaluate and construct vertical coordinate
# Author: Stephanie Westerhuis
# Date: October 2022
#########################################################

# python packages
import numpy as np
import pandas as pd
import iconarray
import xarray as xr
from ipdb import set_trace
import sys
from pathlib import Path
import os
from scipy import spatial


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
        print(f" Given lon: {lon:.3f} vs found lon: {lons[ind]:.3f}")

    return ind


def ind_from_latlon_remapped(lats, lons, lat, lon, verbose=False):
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
    ind_lat = np.argmin(np.abs(lats - lat))
    ind_lon = np.argmin(np.abs(lons - lon))

    if verbose:
        print(f"Closest ind: {ind_lat}, {ind_lon}")
        print(f" Given lat: {lat:.3f} vs found lat: {lats[ind_lat]:.3f}")
        print(f" Given lon: {lon:.3f} vs found lon: {lons[ind_lon]:.3f}")

    return ind_lat, ind_lon


def ind_from_latlon_regular(lats, lons, lat, lon, verbose=False):
    """Find the nearest indices to given location on a regular grid.

    Args:
        lats (2d array):            Latitude grid
        lons (2d array):            Longitude grid
        lat (float):                Latitude of location
        lon (float):                Longitude of location
        verbose (bool, optional):   Print information. Defaults to False.

    Returns:
        ind_lat, ind_lon     Indices of nearest grid point.
    """
    lon_lat = np.array(list(zip(lons.ravel(), lats.ravel())))
    kdtree = spatial.KDTree(lon_lat)
    x_shp = lons.shape[0]
    y_shp = lons.shape[1]
    data_index = np.arange(0, x_shp * y_shp).reshape(x_shp, y_shp)
    distance, i = kdtree.query((lon, lat))
    ind_lat, ind_lon = np.argwhere(data_index == i)[0]

    if verbose:
        print(f"Closest indices in y- and x-direction: {ind_lat}, {ind_lon}")
        print(f" Given lat: {lat:.3f} vs found lat: {lats[ind_lat, ind_lon]:.3f}")
        print(f" Given lon: {lon:.3f} vs found lon: {lons[ind_lat, ind_lon]:.3f}")
        print(f" Distance: {distance:.3f}")

    return ind_lat, ind_lon


def get_poi(loc, lats=None, lons=None, model="icon"):
    """Points of interest for analysis.

    Args:
        lats
        lons
        loc

    Returns:
        pd dataframe

    """

    print("--- Retrieving points of interest")
    all_poi = pd.DataFrame(
        columns=[
            "mtblanc",
            "zrh",
            "pay",
            "visp",
            "ulr",
            "sav",
            "duf",
            "cic",
            "loc",
            "ste",
            "lau",
            "cat",
            "ruc",
        ],
        index=["long_name", "ind", "h_real", "lat", "lon", "left_to_right"],
    )

    all_poi["mtblanc"].long_name = "Mt Blanc"
    all_poi["zrh"].long_name = "Zürich"
    all_poi["pay"].long_name = "Payerne"
    all_poi["visp"].long_name = "Visp"
    all_poi["ulr"].long_name = "Ulrichen"
    all_poi["sav"].long_name = "Savona"
    all_poi["duf"].long_name = "Dufourspitze"
    all_poi["cic"].long_name = "Cicognola"
    all_poi["loc"].long_name = "Locarno"
    all_poi["ste"].long_name = "Steffisburg"
    all_poi["lau"].long_name = "Lausanne"
    all_poi["cat"].long_name = "Catogne"
    all_poi["ruc"].long_name = "Ruchen"

    all_poi["mtblanc"].lat = 45.83267
    all_poi["zrh"].lat = 47.46218
    all_poi["pay"].lat = 46.81291
    all_poi["visp"].lat = 46.29861
    all_poi["ulr"].lat = 46.50568
    all_poi["sav"].lat = 44.27691
    all_poi["duf"].lat = 45.93692
    all_poi["cic"].lat = 45.72350
    all_poi["loc"].lat = 46.16961
    all_poi["ste"].lat = 46.77884
    all_poi["lau"].lat = 46.53447
    all_poi["cat"].lat = 46.0699
    all_poi["ruc"].lat = 47.01078

    all_poi["mtblanc"].lon = 6.86437
    all_poi["zrh"].lon = 8.54458
    all_poi["pay"].lon = 6.94418
    all_poi["visp"].lon = 7.88004
    all_poi["ulr"].lon = 8.30610
    all_poi["sav"].lon = 8.54675
    all_poi["duf"].lon = 7.86675
    all_poi["cic"].lon = 8.61444
    all_poi["loc"].lon = 8.77102
    all_poi["ste"].lon = 7.63525
    all_poi["lau"].lon = 6.58822
    all_poi["cat"].lon = 7.1317
    all_poi["ruc"].lon = 9.00086

    all_poi["mtblanc"].h_real = 4808.0
    all_poi["zrh"].h_real = 422.0
    all_poi["pay"].h_real = 491.0
    all_poi["visp"].h_real = 646.0
    all_poi["ulr"].h_real = 1345.0
    all_poi["sav"].h_real = 0.0
    all_poi["duf"].h_real = 4634.0
    all_poi["cic"].h_real = 197.0
    all_poi["loc"].h_real = 202.0
    all_poi["ste"].h_real = 586.0
    all_poi["lau"].h_real = 415.0
    all_poi["cat"].h_real = 1160.0
    all_poi["ruc"].h_real = 2855.0

    all_poi["mtblanc"].left_to_right = False
    all_poi["zrh"].left_to_right = None
    all_poi["pay"].left_to_right = True
    all_poi["visp"].left_to_right = True
    all_poi["ulr"].left_to_right = False
    all_poi["sav"].left_to_right = None
    all_poi["duf"].left_to_right = False
    all_poi["cic"].left_to_right = True
    all_poi["loc"].left_to_right = None
    all_poi["ste"].left_to_right = False
    all_poi["lau"].left_to_right = False
    all_poi["cat"].left_to_right = False
    all_poi["ruc"].left_to_right = False

    if loc[0] == "all":
        poi = all_poi
    else:
        poi = all_poi[list(loc)]

    # indices of specific locations
    if isinstance(lats, np.ndarray) and isinstance(lons, np.ndarray):
        if "icon_regular" in model:
            retrieve_ind = ind_from_latlon_remapped
        elif "icon" in model:
            retrieve_ind = ind_from_latlon
        elif "cosmo" in model:
            retrieve_ind = ind_from_latlon_regular

        for name, col in poi.items():
            print(f"---  {name}")
            col.ind = retrieve_ind(lats, lons, col.lat, col.lon)

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


def parse_out_dir(out_dir):
    """Parse out_dir input.

    Return figures folders in repo if out_dir = default

    """
    if out_dir == "vc_scratch":
        user = os.getlogin()
        return Path("/scratch/e1000/meteoswiss/scratch", user, "vert_coord/figures")
    else:
        return out_dir


def indices_transect(
    ind, neighbors, left_to_right=True, n_cells_se=35, n_cells_nw=25, verify=False
):

    # triangle collection index
    if left_to_right:
        index_pattern_se = [2, 1] * n_cells_se
    else:
        index_pattern_se = [1, 2] * n_cells_se
    n_all_cells_se = 2 * n_cells_se
    index_list_se = np.empty(n_all_cells_se, dtype=int)

    if left_to_right:
        index_pattern_nw = [1, 2] * n_cells_nw
    else:
        index_pattern_nw = [2, 1] * n_cells_nw
    n_all_cells_nw = 2 * n_cells_nw
    index_list_nw = np.empty(n_all_cells_nw, dtype=int)

    # extend in south-east direction
    new_index = ind
    for e, i in enumerate(index_pattern_se):
        neighbors_new_ind = neighbors[:, new_index]
        new_index = neighbors_new_ind[i] - 1
        index_list_se[e] = new_index

    # extend in north-west direction
    new_index = ind
    for e, i in enumerate(index_pattern_nw):
        neighbors_new_ind = neighbors[:, new_index]
        new_index = neighbors_new_ind[i] - 1
        index_list_nw[e] = new_index
    # reverse order of those indices
    index_list_nw = index_list_nw[::-1]

    # concatenate to one line
    ind_line = np.append(np.append(index_list_nw, np.array(ind)), index_list_se)

    # numbering of indices with respect to origin
    ind_wrt_origin = np.arange(-n_all_cells_nw, n_all_cells_se + 1)

    if verify:
        for e, i in enumerate(ind_line):
            print(f"---- Index at:                      {i}")
            if e % 2 == 0:
                position_next_ind = 2
            else:
                position_next_ind = 1
            print(
                f"Next neighbour index at position {position_next_ind}: {neighbors[position_next_ind,i] - 1 }"
            )

    return ind_line, ind_wrt_origin


def retrieve_lats_lons_hhl_icon(ds):

    try:
        hhl = ds.HEIGHT.values
        lats = ds.clat_1.values
        lons = ds.clon_1.values
    except AttributeError:
        try:
            hhl = ds.HEIGHT.values
            lats = ds.clat.values
            lons = ds.clon.values
        except AttributeError:
            try:
                hhl = ds.HHL.values
                lats = ds.clat.values
                lons = ds.clon.values
            except AttributeError:
                try:
                    hhl = np.array([ds.topography_c.values])
                    lats = ds.clat.values
                    lons = ds.clon.values
                except AttributeError:
                    print(
                        "--- names for 3D height field, latitudes or longitudes unknown!"
                    )
                    sys.exit(0)

    # convert from radians to degrees if necessary
    if max(lats) < 2:
        lats = np.rad2deg(lats)
        lons = np.rad2deg(lons)

    return lats, lons, hhl


def retrieve_vars_print(file_str, model):

    try:
        ds = xr.open_dataset(file_str)
    except FileNotFoundError:
        print(f"!! File does not exist: {file_str}")

    if "icon" in model:
        lats, lons, hhl = retrieve_lats_lons_hhl_icon(ds)

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

    return lats, lons, hhl, hsurf, dz


def retrieve_lats_lons_icon(ds):
    try:
        lats = ds.clat_1.values
        lons = ds.clon_1.values
    except AttributeError:
        try:
            lats = ds.clat.values
            lons = ds.clon.values
        except AttributeError:
            print("Name of latitude and longitude unknown!")
            sys.exit()

    # convert from radians to degrees if necessary
    if max(lats) < 2:
        lats = np.rad2deg(lats)
        lons = np.rad2deg(lons)

    return lats, lons


def retrieve_vars_icon_fcst_u(file):

    try:
        ds = xr.open_dataset(file).squeeze()
    except FileNotFoundError:
        print(f"!! File does not exist: {file_str}")

    lats, lons = retrieve_lats_lons_icon(ds)
    u_wind = ds["U"].values

    return lats, lons, u_wind


def open_icon_ds(file, grid_file):
    ds = iconarray.combine_grid_information(file, grid_file)
    ds_grid = xr.open_dataset(grid_file)

    return ds, ds_grid


def open_ds_regular(file):
    return xr.open_dataset(file)


def retrieve_vars_icon_regular(ds):
    lats = ds.lat.values
    lons = ds.lon.values
    hhl = ds.HSURF.values

    return lats, lons, hhl


def retrieve_vars_cosmo_regular(ds):
    lats = ds.lat_1.values
    lons = ds.lon_1.values
    hsurf = ds.HSURF.values

    return lats, lons, hsurf
