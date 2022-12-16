import iconarray
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ipdb import set_trace
import psyplot.project as psy
import cartopy.crs as ccrs

from utils import indices_transect


def transect_hhl(hhl, neighbors, poi, config, out_dir, n_levels):

    for location in poi:
        loc = poi[location]

        # retrieve indices of cells along 1 straight line
        ind_line, ind_wrt_origin = indices_transect(
            loc.ind, neighbors, n_cells_se=15, n_cells_nw=15
        )

        # create figure
        fig, ax = plt.subplots(1, 1, figsize=(11, 5))

        # transect from A to B
        transect = hhl[:, ind_line]

        # plot surface
        ax.plot(ind_wrt_origin, transect[-1, :], linewidth=2, color="k")

        # plot vertical coordinate surfaces above surface
        l_colors = plt.cm.cividis(np.linspace(0, 1, n_levels))
        for i, level in enumerate(np.arange(2, n_levels + 1)):
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
