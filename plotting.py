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

        transect = hhl[:, ind_line]

        ax.plot(ind_wrt_origin, transect[-1, :], linewidth=2, color="k")

        l_colors = plt.cm.RdPu_r(np.linspace(0, 1, 19))
        for i, level in enumerate(np.arange(2, 20)):
            ax.plot(ind_wrt_origin, transect[-level, :], color=l_colors[i])

        ax.axvline(0, linewidth=0.5, color="grey")

        # plot labelling
        ax.set_title(f"Transect through {loc.long_name}")
        ax.set_xlabel(f"Cells with respect to {loc.long_name}")
        ax.set_ylabel("Altitude [masl]")

        out_name = Path(out_dir, f"transect_{location}_{config}.png")
        plt.savefig(out_name)


def transect_topo(hhl, neighbour_ind, poi, config, out_dir):

    for location in poi:
        loc = poi[location]

        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        # ax.plot(transect[-1, :])

        ## plot labelling
        # ax.set_title(f"Topography starting at {loc.long_name}")
        # ax.set_xlabel("Cells")
        # ax.set_ylabel("Altitude [masl]")

        # out_name = Path(out_dir, f"topography_{location}_{config}.png")
        # plt.savefig(out_name)
