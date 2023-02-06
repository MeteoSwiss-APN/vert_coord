import iconarray
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ipdb import set_trace
import psyplot.project as psy
import cartopy.crs as ccrs
import warnings

warnings.filterwarnings("ignore", message="Shapefile")

from utils import indices_transect
from utils import get_poi
from utils import retrieve_lats_lons_hhl_icon


def transect_hhl(ds, ds_grid, loc, config, out_dir, lev):

    lats, lons, hhl = retrieve_lats_lons_hhl_icon(ds)
    neighbs = ds_grid.neighbor_cell_index.values
    poi = get_poi(loc, lats, lons)

    for location in poi:
        loc = poi[location]

        # retrieve indices of cells along 1 straight line
        ind_line, ind_wrt_origin = indices_transect(
            loc.ind, neighbs, n_cells_se=15, n_cells_nw=15
        )

        # create figure
        fig, ax = plt.subplots(1, 1, figsize=(11, 5))

        # transect from A to B
        transect = hhl[:, ind_line]

        # plot surface
        ax.plot(ind_wrt_origin, transect[-1, :], linewidth=2, color="k")

        # plot vertical coordinate surfaces above surface
        l_colors = plt.cm.cividis(np.linspace(0, 1, lev))
        for i, level in enumerate(np.arange(2, lev + 1)):
            ax.plot(ind_wrt_origin, transect[-level, :], color=l_colors[i])

        # indicate location of poi
        ax.axvline(0, linewidth=0.5, color="grey")

        # plot labelling
        ax.set_title(f"Transect through {loc.long_name}")
        ax.set_xlabel(f"Cells with respect to {loc.long_name}")
        ax.set_ylabel("Altitude [masl]")

        # save
        out_name = Path(out_dir, f"transect_{location}_{config}.png")
        plt.savefig(out_name, dpi=200)
        print(f"Saved as: {out_name}")


def transect_topo(hhl, neighbors, ds, poi, config, out_dir, level=1):

    for location in poi:
        loc = poi[location]
        surf = hhl[-level, :]

        # retrieve indices of cells along 1 straight line
        ind_line, ind_wrt_origin = indices_transect(
            loc.ind, neighbors, n_cells_se=15, n_cells_nw=15
        )

        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        # fig, ax2 = plt.subplots(1, 1, figsize=(10, 5))
        fig = plt.figure(figsize=(7, 2.7), tight_layout=False)
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.2, top=0.95)
        gs = fig.add_gridspec(
            nrows=1,
            ncols=10,
        )
        ax1 = fig.add_subplot(gs[0, 0:3], projection=ccrs.PlateCarree())
        ax2 = fig.add_subplot(gs[0, 4:])

        # 1st part of figure: 2d orography
        ##################################
        mask = np.ones(len(surf))
        mask[ind_line] = np.nan
        ds = ds.assign(HSURF_masked=ds["HSURF"] * mask)
        ds.HSURF_masked.encoding["coordinates"] = "clat clon"
        mapplot = psy.plot.mapplot(
            ds,
            name="HSURF_masked",
            map_extent=[loc.lon - 0.7, loc.lon + 0.7, loc.lat - 0.6, loc.lat + 0.6],
            xgrid=False,
            ygrid=False,
            cmap="viridis",
            ax=ax1,
        )

        # 2nd part of figure: vertical transect
        #######################################

        # transect from A to B
        transect_surf = surf[ind_line]

        # plot surface
        ax2.plot(
            ind_wrt_origin, transect_surf, linewidth=2, color="darkorange", label="orig"
        )

        # add additional line to indicate every 2nd grid cell
        transect_surf_2 = transect_surf[np.arange(0, len(ind_line), 2)]
        ind_wrt_origin_2 = ind_wrt_origin[np.arange(0, len(ind_line), 2)]
        ax2.plot(
            ind_wrt_origin_2,
            transect_surf_2,
            linewidth=0.4,
            color="darkblue",
            label="every 2nd cell",
        )

        # indicate location of poi
        ax2.axvline(0, linewidth=0.5, color="grey")

        # plot labelling
        ax2.set_xlabel(f"Cells with respect to {loc.long_name}")
        ax2.set_ylabel("Altitude [masl]")

        # legend
        ax2.legend()

        # save
        out_name = Path(out_dir, f"topo_{location}_{config}.png")
        plt.savefig(out_name, dpi=200)
        print(f"Saved as: {out_name}")


def mapplot_coord_surf(ds, config, out_dir, lev, loc, radius, vmin, vmax):

    fig = plt.figure(figsize=(12, 7), tight_layout=True)
    ax1 = fig.add_subplot(projection=ccrs.PlateCarree())

    # get level number, find variable name for height field
    try:
        nlev = ds.HHL.values.shape[0]
        name_height = "HHL"
    except AttributeError:
        nlev = ds.HEIGHT.values.shape[0]
        name_height = "HEIGHT"
    zz = nlev - lev

    if loc[0] == "all":
        surf_map = ds.psy.plot.mapplot(
            ax=ax1,
            name=name_height,
            cticksize="small",
            borders=True,
            lakes=True,
            rivers=True,
            projection="robin",
            title=f"Elevation of {lev}. coordinate surface",
            clabel="m asl",
            z=zz,  # specify level
        )
        out_name = Path(out_dir, f"altitude_{lev}_coordinate_surface_{config}.png")
        plt.savefig(out_name)
        print(f"Saved as: {out_name}")

    else:
        poi = get_poi(loc)

        for name, col in poi.items():
            print(f"--- Plotting {lev}. level around {name}...")

            lonmin = col.lon - radius * 0.01
            lonmax = col.lon + radius * 0.01
            latmin = col.lat - radius * 0.01
            latmax = col.lat + radius * 0.01

            surf_map = ds.psy.plot.mapplot(
                ax=ax1,
                name=name_height,
                xgrid=None,
                ygrid=None,
                map_extent=[lonmin, lonmax, latmin, latmax],
                cmap="terrain",
                cticksize="small",
                projection="robin",
                title=f"Elevation of {lev}. coordinate surface around {col.long_name}",
                clabel="m asl",
                z=zz,  # specify specific level
                bounds={"method": "minmax", "vmin": vmin, "vmax": vmax},
            )

            out_name = Path(
                out_dir, f"altitude_{lev}_coordinate_surface_{name}_{config}.png"
            )
            plt.savefig(out_name)
            print(f"Saved as: {out_name}")


def profile_dz(dz, hhl, poi, loc, config, out_dir):
    for location in poi:
        loc = poi[location]

        fig, ax = plt.subplots(1, 1, figsize=(5, 7))

        ax.plot(dz[:, loc.ind], hhl[1:, loc.ind])

        # plot labelling
        ax.set_title(f"delta_z at {loc.long_name}")
        ax.set_xlabel(f"Level thickness [m]")
        ax.set_ylabel(f"Altitude [masl]")

        # save
        out_name = Path(out_dir, f"ddz_{location}_{config}.png")
        plt.tight_layout()
        plt.savefig(out_name, dpi=200)
        print(f"Saved as: {out_name}")


#    """Plot dz vs altitude.
#
#    Args:
#        dz (np.array): level thickness
#        poi (pd.Dataframe): points of interest
#        model (str): model name
#        exp (str): identifier for specific configuration
#        out_dir (str): output directory
#
#    """
#    # plot "nominal" dz, i.e. at location of savona
#    ii = poi["sav"].ind
#    dz_nominal = dz[::-1, ii]
#    n_lev = len(dz_nominal)
#
#    # cut into lower and upper region
#    cut = 22
#    dz_lower = dz_nominal[:cut]
#    dz_upper = dz_nominal[(cut - 1) :]
#
#    ### dz vs level-index ###
#    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 6))
#
#    # plot all levels
#    ax1.plot(np.arange(1, cut + 1), dz_lower, marker="o", color="purple")
#    ax1.plot(np.arange(cut, n_lev + 1), dz_upper, marker="o", color="yellowgreen")
#    ax1.set_xlabel("# Level")
#    ax1.set_ylabel("dz [m]")
#    ax1.set_title("All levels")
#    ax1.set_ylim(0, 1200)
#    ax1.set_xlim(0, n_lev)
#
#    # plot only lower atmosphere
#    ax2.plot(np.arange(1, cut + 1), dz_lower, marker="o", color="purple")
#    ax2.set_xlabel("# Level")
#    ax2.set_title("Lower atmosphere")
#    ax2.set_ylim(min(dz_lower) - 2, max(dz_lower) + 5)
#    ax2.set_xlim(0, cut)
#
#    Path(out_dir).mkdir(parents=True, exist_ok=True)
#    out_name = Path(out_dir, f"dz_vs_level_{model}_{exp}.png")
#    plt.tight_layout()
#    plt.savefig(out_name)
#    print(f"Saved as: {out_name}")
#
#    plt.clf()
#
#    ### altitude vs dz ###
#    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 8))
#
#    # calculate cumulative sum of dz
#    dz_cumsum = np.cumsum(dz_nominal)
#
#    # plot all levels
#    ax1.plot(dz_lower, dz_cumsum[:cut], marker="o", color="purple")
#    ax1.plot(dz_upper, dz_cumsum[(cut - 1) :], marker="o", color="yellowgreen")
#    ax1.set_xlabel("dz [m]")
#    ax1.set_ylabel("Altitude [m asl]")
#    ax1.set_title("All levels")
#    ax1.set_xlim(0, 1200)
#    ax1.set_ylim(0, 20500)
#
#    # plot only lower atmosphere
#    ax2.plot(dz_lower, dz_cumsum[:cut], marker="o", color="purple")
#    ax2.set_xlabel("dz [m]")
#    ax2.set_title("Lower atmosphere")
#    ax2.set_xlim(min(dz_lower) - 2, max(dz_lower) + 2)
#    ax2.set_ylim(0, dz_cumsum[(cut - 1)])
#
#    Path(out_dir).mkdir(parents=True, exist_ok=True)
#    out_name = Path(out_dir, f"altitude_vs_dz_{model}_{exp}.png")
#    plt.tight_layout()
#    plt.savefig(out_name)
#    print(f"Saved as: {out_name}")
#
#


def transect_topo_regular(lats, lons, hsurf, poi, radius, config, out_dir):

    if "reg" in config:
        radius = radius * 10

    for location in poi:
        loc = poi[location]

        # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        # fig, ax2 = plt.subplots(1, 1, figsize=(10, 5))
        fig = plt.figure(figsize=(7, 2.7), tight_layout=False)
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.2, top=0.95)
        gs = fig.add_gridspec(
            nrows=1,
            ncols=10,
        )
        ax1 = fig.add_subplot(gs[0, 0:3])  # , projection=ccrs.PlateCarree())
        ax2 = fig.add_subplot(gs[0, 4:])

        # 1st part of figure: 2d orography
        ##################################
        ind_lat = loc.ind[0]
        ind_lon = loc.ind[1]
        lonmin = ind_lon - radius
        lonmax = ind_lon + radius
        latmin = ind_lat - radius
        latmax = ind_lat + radius
        regional_hsurf = hsurf[latmin:latmax, lonmin:lonmax]
        vmin = regional_hsurf.min()
        vmax = regional_hsurf.max()
        ax1.imshow(regional_hsurf, origin="lower", vmin=vmin, vmax=vmax, cmap="terrain")
        ax1.scatter(radius, radius, marker="s", color="red", s=15)

        # 2nd part of figure: vertical transect
        #######################################

        # plot surface
        ax2.plot(np.arange(-radius, radius), hsurf[ind_lat, lonmin:lonmax])
        # ax2.plot(
        #    ind_wrt_origin, transect_surf, linewidth=2, color="darkorange", label="orig"
        # )

        ## indicate location of poi
        ax2.axvline(0, linewidth=0.5, color="grey")

        ## plot labelling
        ax2.set_xlabel(f"Cells to the west and east of {loc.long_name}")
        ax2.set_ylabel("Altitude [masl]")

        ## legend
        # ax2.legend()

        # save
        out_name = Path(out_dir, f"topo_regular_{location}_{config}.png")
        plt.savefig(out_name, dpi=200)
        print(f"Saved as: {out_name}")
