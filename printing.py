import numpy as np
import iconarray
from ipdb import set_trace

from utils import get_min_max
from utils import n_sum_up_to


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
