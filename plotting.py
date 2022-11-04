import iconarray
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ipdb import set_trace


def lineplot_transect(hhl, neighbour_ind, poi, config, out_dir):

    for location in poi:
        loc = poi[location]

        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        # triangle collection index
        index_pattern_se = [2, 1] * 35
        index_list_se = np.empty(70, dtype=int)
        index_pattern_nw = [1, 2] * 15
        index_list_nw = np.empty(30, dtype=int)

        # extend in south-east direction
        new_index = loc.ind
        for e, i in enumerate(index_pattern_se):
            neighbors = neighbour_ind[:, new_index]
            new_index = neighbors[i] - 1
            index_list_se[e] = new_index

        # extend in northwest direction
        new_index = loc.ind
        for e, i in enumerate(index_pattern_nw):
            neighbors = neighbour_ind[:, new_index]
            new_index = neighbors[i] - 1
            index_list_nw[e] = new_index

        ind_line = np.append(np.append(index_list_nw, np.array(loc.ind)), index_list_se)

        transect = hhl[:, ind_line]

        for i in range(1, 20):
            ax.plot(transect[-i, :])

        # plot labelling
        ax.set_title(f"Transect through {loc.long_name}")
        ax.set_xlabel("Cells")
        ax.set_ylabel("Altitude [masl]")

        out_name = Path(out_dir, f"transect_{location}_{config}.png")
        plt.savefig(out_name)
