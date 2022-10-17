# Load modules
import cmcrameri.cm as cmc
import cartopy.feature as cf                                                                                                        
from pathlib import Path
import psyplot.project as psy
import numpy as np
import click
import sys

from iconarray.plot import formatoptions # import plotting formatoptions (for use with psyplot)
import iconarray as iconvis # import self-written modules from iconarray



# Define functions used
'''        
    Function fill_geofac_div, PSEUDO CODE:     
    LOOP over triangles (T)
        LOOP over edges (E)
            COEFF(E) = PRIMAL_EDGE_LENGTH(E)*EDGE_ORIENTATION(E)/AREA(T)
    
'''
def fill_geofac_div(grid, Ntriangles = 100):
    # assume that cell_type = 3
    cell_type = 3
    geofac_div = np.zeros((Ntriangles, cell_type))
    for it in range(Ntriangles):
        for ie in range(cell_type):
            geofac_div[it,ie] = grid.edge_length.values[grid.edge_of_cell[ie,it]] * \
            grid.edge_system_orientation.values[grid.vertex_of_cell[ie,it]]/        \
            grid.cell_area.values[it]
    return geofac_div

''' 
    Function fill_geofac_n2s, PSEUDO CODE:
    LOOP over triangles (T)
        LOOP over edges (E)
            T1, T2 <- two triangles that share E
            IF T = T1
                COEFF(1) = COEFF(1) - DIV(E)/DUAL_EDGE_LENGTH(E)
            ELIF T = T2
                COEFF(1) = COEFF(1) + DIV(E)/DUAL_EDGE_LENGTH(E)
            LOOP over edges (E2)
                TE <- neighbor (T, E2) ! neighbor of T that shares E2
                IF TE = T1
                    COEFF(E2 + 1) = COEFF(E2 + 1) - DIV(E)/DUAL_EDGE_LENGTH(E)
                ELIF TE = T2
                    COEFF(E2 + 1) = COEFF(E2 + 1) + DIV(E)/DUAL_EDGE_LENGTH(E)
'''
def fill_geofac_n2s(grid, geofac_div, Ntriangles = 100):
    # assume that cell_type = 3
    cell_type = 3
    geofac_n2s = np.zeros((Ntriangles, cell_type + 1))
    for it in range(Ntriangles):
        for ie in range(cell_type):
            [it1, it2] = grid.adjacent_cell_of_edge.values[:,ie]
            if it == it1:
                geofac_n2s[it,0] = geofac_n2s[it,0] - geofac_div[it,ie]/grid.dual_edge_length.values[ie]
            elif it == it2:
                geofac_n2s[it,0] = geofac_n2s[it,0] + geofac_div[it,ie]/grid.dual_edge_length.values[ie]
            for ie2 in range(cell_type):
                ites = grid.adjacent_cell_of_edge.values[:, ie2] # triangles that share E2
                itns = grid.neighbor_cell_index.values[:,it]     # triangles that are neighbor of T
                ite_ar = np.intersect1d(itns[itns != -1], ites[ites != -1], assume_unique = True)
                if ite_ar:
                    ite = ite_ar[0]
                    if ite == it1:
                        geofac_n2s[it,ie2+1] = geofac_n2s[it,ie2+1] - geofac_div[it,ie]/grid.dual_edge_length.values[ie]
                    elif ite == it2:
                        geofac_n2s[it,ie2+1] = geofac_n2s[it,ie2+1] + geofac_div[it,ie]/grid.dual_edge_length.values[ie]           
    return geofac_n2s
'''
    Function nabla2_scalar
    Returns the laplacian of a field that is defined over all cell points.
    Inputs:  grid         - the grid data
             psi_c        - the input field
             geofac_n2s   - the geometrical factors applied on the stencil, 
                                            output of fill_geofac_n2s
             Ntriangles   - the number of cells on which the computation is done.
    Outputs: nabla2_psi_c - the laplacian of psi_c 
'''
def nabla2_scalar(grid, psi_c, geofac_n2s, Ntriangles):
    nabla2_psi_c = np.zeros((Ntriangles))
    for it in range(Ntriangles):
        nabla2_psi_c[it] = psi_c[it]*geofac_n2s[it,0] +            \
        psi_c[grid.neighbor_cell_index[0, it]]*geofac_n2s[it,1] +  \
        psi_c[grid.neighbor_cell_index[1, it]]*geofac_n2s[it,2] +  \
        psi_c[grid.neighbor_cell_index[2, it]]*geofac_n2s[it,3]
    return nabla2_psi_c

'''
    Taken from the ICON code in atm_dyn_iconam/mo_nh_init_utils.f90
'''
def compute_smooth_topo(grid, topo, geofac_n2s, niter = 25, Ntriangles = 100):
    smooth_topo = np.zeros_like(topo)
    smooth_topo[:Ntriangles] = topo[:Ntriangles]
    for i in range(niter):
        nabla2_topo = nabla2_scalar(grid, topo, geofac_n2s, Ntriangles)
        smooth_topo[:Ntriangles] = smooth_topo[:Ntriangles] + \
        0.125 * nabla2_topo[:Ntriangles]*grid.cell_area.values[:Ntriangles]
        print('Finished iteration : {}'.format(i))
    return smooth_topo
    

@click.command()
@click.option(
    "--file",
    default="/store/mch/msopr/swester/teamx/tdf_2019091212/output/19091212/lfff00000000c.nc",
    help="Netcdf file containing HSURF to be smoothed.",
    type=str,
)
@click.option(
    "--grid_file",
    default="/store/mch/msopr/swester/grids/alps_R19B08/alps_DOM01.nc",
    help="Netcdf file containing grid information.",
    type=str,
)
@click.option(
    "--div_file",
    default="",
    help="File where the divergence coefs. are stored. By default they are computed in this script.",
    type=str,
)
@click.option(
    "--n2s_file",
    default="",
    help="File where the n2s coefs. are stored. By default they are computed in this script.",
    type=str,
)
@click.option(
    "--out_dir",
    default="/scratch/snx3000/gvanpary/vert_coord/figures",
    type=str,
)
@click.option(
    "--save_div",
    is_flag=True,
    default=False,
)
@click.option(
    "--save_n2s",
    is_flag=True,
    default=False,
)
@click.option(
    "--Ntriangles",
    type=int,
    default=1000,
    help="Number of triangles on which to do the computation, from index 0 to Ntriangles."
)
@click.option(
    "--Niter",
    type=int,
    default=25,
    help="Number of smoothing iterations."
)
def smooth_topo(
    file,
    grid_file,
    Ntriangles,
    Niter,
    save_div,
    save_n2s,
    div_file,
    n2s_file,
    out_dir,
):
    # Need to add the option to load div and n2s from file!
    print('Smoothing topology from {file}')
    # Load data
    grid = psy.open_dataset(grid_file)
    data = psy.open_dataset(file)
    if div_file == "":
        div = fill_geofac_div(grid, Ntriangles)
        print('Finished computing divergence coefficients')
    else: 
        #div =
        # implement loading div file (pickle ?)
        print('Loaded divergence coefficients')
    
    if n2s_file == "": 
        n2s = fill_geofac_n2s(grid, div, Ntriangles)
        print('Finished computing n2s coefficients')
    else:
        #div =
        # implement loading div file (pickle ?)
        print('Loaded n2s coefficients')
    smooth_HSURF = compute_smooth_topo(grid, data.HSURF, n2s, Niter, Ntriangles)
    print('Finished smoothing topology')
    # implement saving (using pickle ?)
if __name__ == "__main__":
    smooth_topo()