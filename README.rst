==========
vert_coord
==========

This repo gathers a variety of python scripts to create and evaluate a vertical coordinate for ICON.

.. contents:: **Content**

------------
Installation
------------

This package does not yet comply with high-standard CICD.

Simply init the dedicated conda env with:

``conda env create -f environment.yml``

and later update with:

``conda env update --file environment.yml``

------------------------
Usage of evaluate_hhl.py
------------------------
- ``lev``: Indicate number of levels to be printed
- ``loc``: Locations
- ``grid_file``: ICON grid file
- ``file``: constants file containing HHL
- ``model``: model, either icon or cosmo

Print altitude of levels at specific locations
----------------------------------------------

``python evaluate_hhl.py --print_hhl --lev 35 --loc mtblanc --loc sav --loc zrh --loc ulr --grid_file /store/s83/swester/vert_coord_files/icon-1-alps/alps_DOM01.nc --file /store/s83/swester/vert_coord_files/icon-1-alps/const_sleve.nc``

Print maximum elevation difference between adjacent cells
---------------------------------------------------------
``python evaluate_hhl.py --print_max_dzdc --lev 1 --grid_file /store/s83/swester/vert_coord_files/icon-1-alps/alps_DOM01.nc --file /store/s83/swester/vert_coord_files/icon-1-alps/const_sleve.nc``

-------------------------
Usage of construct_hhl.py
-------------------------
Constructs HHL field from an indicated surface elevation field and other settings. Currently, the script only calculates vct_a and vct_b.

``python construct_hhl.py --file /store/s83/swester/vert_coord_files/icon-1-alps/external_parameter_icon_alps_R19B08_mch.nc --n_levels 81 --h_flat 16000 --top_height 22000 --stretch_fac 0.65 --type_vct_a 2nd_order --type_vct_b linear``


---------------------------------------------------------
plot vertical slices from A to B with losvec_slices.ipynb
---------------------------------------------------------





