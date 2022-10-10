#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%%
z_cosmo = np.array(
    [
        22000.0000,
        21187.5195,
        20396.2500,
        19625.8906,
        18876.1289,
        18146.6797,
        17437.2188,
        16747.4492,
        16077.0781,
        15425.7969,
        14793.2969,
        14179.2891,
        13583.4688,
        13005.5273,
        12445.1875,
        11902.1289,
        11376.0586,
        10866.6875,
        10373.7070,
        9896.8398,
        9435.7695,
        8990.2070,
        8559.8672,
        8144.4492,
        7743.6602,
        7357.2070,
        6984.7969,
        6626.1484,
        6280.9688,
        5948.9570,
        5629.8281,
        5323.3086,
        5029.0898,
        4746.9102,
        4476.4570,
        4217.4570,
        3969.6399,
        3732.7000,
        3506.3701,
        3290.3601,
        3084.3999,
        2888.2000,
        2701.4800,
        2523.9600,
        2355.3799,
        2195.4500,
        2043.8999,
        1900.4500,
        1764.8398,
        1636.7800,
        1516.0200,
        1402.2700,
        1295.2800,
        1194.7700,
        1100.4800,
        1012.1399,
        929.4900,
        852.2800,
        780.2300,
        713.0898,
        650.6099,
        592.5200,
        538.5698,
        488.5198,
        442.1099,
        399.0898,
        359.2100,
        322.2400,
        287.9299,
        256.0398,
        226.3300,
        198.5700,
        172.5400,
        148.0000,
        124.7300,
        102.5100,
        81.1300,
        60.3900,
        40.0700,
        20.0000,
        0.0000,
    ]
)[::-1]

#%%
dz_cosmo = z_cosmo[1:] - z_cosmo[:-1]

#%%
levels_z = np.arange(1, len(z_cosmo) + 1)
levels_dz = np.arange(1, len(dz_cosmo) + 1)

# levels_z = np.arange(len(z_cosmo))
# levels_dz = np.arange(len(dz_cosmo))

#%%
fix, ax = plt.subplots(figsize=(4, 5))
ax.scatter(levels_dz, dz_cosmo, marker="o")
ax.set_yscale("log")

# %%
def expo(x, a, b):
    return a * x**b


# %%
pars, cov = curve_fit(
    f=expo, xdata=levels, ydata=dz_cosmo, p0=[0, 2], bounds=(-100, 100)
)
print(f"a = {pars[0]:.2f}, b = {pars[1]:.2f}")

# %%

fix, ax = plt.subplots(figsize=(6, 8))
ax.scatter(levels, dz_cosmo, marker="o", label="orig")
ax.set_yscale("log")
ax.plot(levels, expo(levels, *pars), linestyle="--", color="r", label="fitted")
plt.legend()
ax.set_xlabel("Level")
ax.set_ylabel("dz")

# %%
for i in range(80):
    print(f"{dz_cosmo[i]:.2f} vs {expo(levels, *pars)[i]:.2f}")
