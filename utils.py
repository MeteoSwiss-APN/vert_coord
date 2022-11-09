# utilities to evaluate and construct vertical coordinate
# Author: Stephanie Westerhuis
# Date: October 2022
#########################################################

# python packages
import numpy as np
import pandas as pd
import iconarray
from ipdb import set_trace
import sys
from pathlib import Path
import os


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
        columns=["mtblanc", "zrh", "pay", "visp", "ulr", "sav", "duf", "cic", "ste"],
        index=["long_name", "ind", "h_real", "lat", "lon"],
    )

    poi["mtblanc"].long_name = "Mt Blanc"
    poi["zrh"].long_name = "Zürich"
    poi["pay"].long_name = "Payerne"
    poi["visp"].long_name = "Visp"
    poi["ulr"].long_name = "Ulrichen"
    poi["sav"].long_name = "Savona"
    poi["duf"].long_name = "Dufourspitze"
    poi["cic"].long_name = "Cicognola"
    poi["ste"].long_name = "Steffisburg"

    poi["mtblanc"].lat = 45.83267
    poi["zrh"].lat = 47.46218
    poi["pay"].lat = 46.81291
    poi["visp"].lat = 46.29861
    poi["ulr"].lat = 46.50568
    poi["sav"].lat = 44.27691
    poi["duf"].lat = 45.93692
    poi["cic"].lat = 45.72350
    poi["ste"].lat = 46.77884

    poi["mtblanc"].lon = 6.86437
    poi["zrh"].lon = 8.54458
    poi["pay"].lon = 6.94418
    poi["visp"].lon = 7.88004
    poi["ulr"].lon = 8.30610
    poi["sav"].lon = 8.54675
    poi["duf"].lon = 7.86675
    poi["cic"].lon = 8.61444
    poi["ste"].lon = 7.63525

    # indices of specific locations
    ind_mtblanc = ind_from_latlon(lats, lons, 45.83267, 6.86437)
    ind_zrh = ind_from_latlon(lats, lons, 47.46218, 8.54458)
    ind_pay = ind_from_latlon(lats, lons, 46.81291, 6.94418)
    ind_visp = ind_from_latlon(lats, lons, 46.29861, 7.88004)
    ind_ulr = ind_from_latlon(lats, lons, 46.50568, 8.30610)
    ind_sav = ind_from_latlon(lats, lons, 44.27691, 8.54675)
    ind_duf = ind_from_latlon(lats, lons, 45.93692, 7.86675)
    ind_cic = ind_from_latlon(lats, lons, 45.72350, 8.61444)
    ind_ste = ind_from_latlon(lats, lons, 46.77884, 7.63525)

    poi["mtblanc"].ind = ind_mtblanc
    poi["zrh"].ind = ind_zrh
    poi["pay"].ind = ind_pay
    poi["visp"].ind = ind_visp
    poi["ulr"].ind = ind_ulr
    poi["sav"].ind = ind_sav
    poi["duf"].ind = ind_duf
    poi["cic"].ind = ind_cic
    poi["ste"].ind = ind_ste

    poi["mtblanc"].h_real = 4808.0
    poi["zrh"].h_real = 422.0
    poi["pay"].h_real = 491.0
    poi["visp"].h_real = 646.0
    poi["ulr"].h_real = 1345.0
    poi["sav"].h_real = 0.0
    poi["duf"].h_real = 4634.0
    poi["cic"].h_real = 197.0
    poi["ste"].h_real = 586.0

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


def indices_transect(ind, neighbors, n_cells_se=35, n_cells_nw=25, verify=False):

    # triangle collection index
    index_pattern_se = [2, 1] * n_cells_se
    n_all_cells_se = 2 * n_cells_se
    index_list_se = np.empty(n_all_cells_se, dtype=int)

    index_pattern_nw = [1, 2] * n_cells_nw
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
