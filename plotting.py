import iconarray
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ipdb import set_trace

from utils import indices_transect


def transect_hhl(hhl, neighbors, poi, config, out_dir):

    for location in poi:
        loc = poi[location]

        # retrieve indices of cells along 1 straight line
        ind_line, ind_wrt_origin = indices_transect(loc.ind, neighbors)

        # create figure
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        # transect from A to B
        transect = hhl[:, ind_line]

        # plot surface
        ax.plot(ind_wrt_origin, transect[-1, :], linewidth=2, color="k")

        # plot vertical coordinate surfaces above surface
        l_colors = plt.cm.RdPu_r(np.linspace(0, 1, 19))
        for i, level in enumerate(np.arange(2, 20)):
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


def transect_topo(hhl, neighbors, poi, config, out_dir):

    for location in poi:
        loc = poi[location]

        # retrieve indices of cells along 1 straight line
        ind_line, ind_wrt_origin = indices_transect(loc.ind, neighbors)

        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        # transect from A to B
        transect_surf = hhl[-1, ind_line]

        # plot surface
        ax.plot(ind_wrt_origin, transect_surf, linewidth=2, color="darkorange")

        # indicate location of poi
        ax.axvline(0, linewidth=0.5, color="grey")

        # plot labelling
        ax.set_title(f"Transect of topography through {loc.long_name}")
        ax.set_xlabel(f"Cells with respect to {loc.long_name}")
        ax.set_ylabel("Altitude [masl]")

        # save
        out_name = Path(out_dir, f"topo_{location}_{config}.png")
        plt.savefig(out_name, dpi=200)
        print(f"Saved as: {out_name}")
