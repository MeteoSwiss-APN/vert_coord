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

Print altitude of levels at specific locations
----------------------------------------------
- lev: Indicate number of levels to be printed
- loc: Locations
- grid_file: ICON grid file
- file: constants file containing HHL

``python evaluate_hhl.py --print_hhl --lev 35 --loc mtblanc --loc sav --loc zrh --loc ulr --grid_file /store/s83/swester/vert_coord_files/icon-1-alps/alps_DOM01.nc --file /store/s83/swester/vert_coord_files/icon-1-alps/const_sleve.nc``



-------------------------
Usage of construct_hhl.py
-------------------------

------------------------------------------------------------------------
plot vertical slices from A to B with losvec_slices.ipynb
------------------------------------------------------------------------





