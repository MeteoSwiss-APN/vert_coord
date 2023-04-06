# Illustrate difference between real orography (ASTER 30m resolved data)
#  and a 1km NWP orography (COSMO1)

# python packages
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from ipdb import set_trace
from pathlib import Path

# homemade packages
from utils import ind_from_latlon

# files
f_aster = "/store/c2sm/extpar_raw_data/topo/aster/ASTER_orig_T031.nc"
f_cosmo = "/store/s83/swester/grids/const_modinterim.nc"
out_dir = "/scratch/e1000/meteoswiss/scratch/swester/vert_coord/oro_aster_model"

# domains
#latmin = 46.4
#latmax = 46.8
#lonmin = 7.7
#lonmax = 8.2
latmin = 46.45
latmax = 46.75
lonmin = 7.9
lonmax = 8.3

levels = np.arange(400, 4200, 200)
d_aster = xr.open_dataset(f_aster).squeeze()
d_cosmo = xr.open_dataset(f_cosmo).squeeze()

lats_aster = d_aster.lat
lons_aster = d_aster.lon
oro_aster = d_aster.Z

lats_cosmo = d_cosmo.lat_1
lons_cosmo = d_cosmo.lon_1
oro_cosmo = d_cosmo.HSURF

# plotting
##########

fig, (ax_a, ax_c) = plt.subplots(ncols=2, figsize=(10, 5.5), sharey=True)

# aster plot
lons_aster_zoom = lons_aster.sel(lon=slice(lonmin, lonmax))
lats_aster_zoom = lats_aster.sel(lat=slice(latmax, latmin))
oro_aster_zoom = oro_aster.sel(lat=slice(latmax, latmin), lon=slice(lonmin, lonmax))
aster = ax_a.contourf(
    lons_aster_zoom,
    lats_aster_zoom,
    oro_aster_zoom,
    cmap="terrain",
    levels=levels,
    #algorithm='serial'
)

ax_a.plot(7.85, 46.68, "o", color="r", markersize=8)
ax_a.text(7.86, 46.67, "Interlaken")
ax_a.plot(7.96, 46.54, "o", color="orangered", markersize=8)

ax_a.set_xlim(lonmin, lonmax)
ax_a.set_ylim(latmin, latmax)

# cosmo plot
cosmo = ax_c.contourf(
    lons_cosmo,
    lats_cosmo,
    oro_cosmo,
    cmap="terrain",
    levels=levels,
)

ax_c.plot(7.85, 46.68, "o", color="r", markersize=8)
ax_c.plot(7.96, 46.54, "o", color="orangered", markersize=8)
ax_c.text(7.97, 46.53, "Jungfrau")

ax_c.set_xlim(lonmin, lonmax)
ax_c.set_ylim(latmin, latmax)

fig.subplots_adjust(
    wspace=0.06, hspace=0.06, left=0.05, right=0.98, top=0.98, bottom=0.2
)
cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.05])
cbar = fig.colorbar(cosmo, cax=cbar_ax, orientation="horizontal", shrink=0.5)

out_name = Path(out_dir, f"oro_aster_vs_cosmo.png")
plt.savefig(out_name, dpi=200)
print(f"Saved as: {out_name}")
plt.clf()
