## Compare aerosol size distributions

```python
>>> from __future__ import print_function, division
...
>>> import numpy as np
>>> import pandas as pd
>>> import xray
...
>>> %matplotlib inline
>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> sns.set(style='ticks', context="talk")
```

Read in a dataset which has been processed with the [`calc_aerosol.py`](calc_aerosol.py) script to extract size distribution parameters.

```python
>>> data = xray.open_dataset("new_aerosol_dists.new_aerosols.nc",
...                          decode_times=False)
...
>>> # Nudge times to the year 2000
... times = data.coords.to_dataset().time
>>> times += 2000.*365
>>> data = xray.decode_cf(data)
...
>>> # Global troposphere slice for quick ref
... global_tropo = data.sel(lev=slice(500, 1100), lat=slice(-80, 80))
>>> global_tropo = global_tropo.isel(time=-1)
```

Quick access to all combinations of data slices

```python
>>> from itertools import product
>>> from plot_funcs import lats, levels
...
>>> def all_subsets(supersets=True):
...     for lev, lev_slice in levels.iteritems():
...         if supersets:
...             yield ((lev, lev_slice), ("global", slice(-80, 80)))
...
...         for lat, lat_slice in lats.iteritems():
...             if lat == 'global': continue
...             yield (lev, lev_slice), (lat, lat_slice)
```

Try plotting distributions just for the stratosphere data, using the specified time subset.

```python
>>> from plot_funcs import all_modes, aerosol_boxplots
>>> from marc_analysis import save_figure
...
>>> plt.ioff()
...
>>> for subset in all_subsets():
...     (lev, lev_slice), (lat, lat_slice) = subset
...     print(lev, lat)
...
...     data_subset = data.sel(lev=lev_slice, lat=lat_slice)
...     data_subset = data_subset.isel(time=-1)
...
...     id_str = "{}_{}".format(lev, lat)
...
...     for mode in all_modes:
...         print(mode)
...         if mode != 'MBS': continue
...         fig = aerosol_boxplots(mode, data_subset)
...         save_figure("{}_{}".format(mode, id_str),
...                     fig=fig, qual='quick')
...         plt.close(fig)
```

Percentile data tables

```python
>>> from plot_funcs import mode_percentiles
...
>>> n_pcts = mode_percentiles(global_tropo, 'n')
>>> mu_pcts = mode_percentiles(global_tropo, 'mu')
```

```python
>>> print(n_pcts.swapaxes(0, 1)[['5', '95']])
>>> print(mu_pcts.swapaxes(0, 1)[['5', '95']])
```

What percentile does a given number concentration fall into for each aerosol mode?

```python
>>> from scipy.stats import percentileofscore
>>> from plot_funcs import all_modes
...
>>> n = 0.1
>>> data_set = global_tropo
...
>>> print("% of size distributions where N < {} for each mode".format(n))
>>> for mode in all_modes:
...     print("{:10s}".format(mode), end=" - ")
...     aer_data = data_set['n' + mode].values.ravel()
...
...     pctile = percentileofscore(aer_data, n)
...     print(pctile)
```

Compile a faceted plot showing boxplots for each of the different spatial distributions for each mode, separately

```python
>>> %autoreload
>>> import plot_funcs
>>> cb = plot_funcs.compare_boxplots
...
>>> from marc_analysis import save_figure
...
>>> plt.ioff()
>>> for (mode, param) in product(plot_funcs.all_modes, ['mu', 'n']):
...     if ((mode.startswith("DST") or mode.startswith("SSLT"))
...         and (param == 'mu')): continue
...     g = cb(mode,param,all_subsets,data)
...     save_figure("{}_{}_comparison".format(mode, param),
...                 fig=g.fig, qual='quick')
...     plt.close(g.fig)
```

Diffrences between modes for low-level (near-sfc) aerosol in continental and maritime regimes, at all latitude combinations.

```python
>>> from marc_analysis import save_figure
>>> from marc_analysis.analysis import extract_feature
>>> import plot_funcs
>>> cmc = plot_funcs.compare_maritime_vs_continent
...
>>> from itertools import product
>>> sns.set(style='ticks', context="talk")
...
>>> ocean_mask = extract_feature(data.isel(time=0, lev=0)['T'])
>>> data['ocean'] = (('lat', 'lon', ), ocean_mask)
...
>>> plt.ioff()
>>> for (mode, param) in product(plot_funcs.all_modes, ['mu', 'n']):
...     if ((mode.startswith("DST") or mode.startswith("SSLT"))
...         and (param == 'mu')): continue
...     g = cmc(mode,param,data)
...     save_figure("{}_{}_ocean_vs_land_comparison".format(mode, param),
...                 fig=g.fig, qual='quick')
...     plt.close(g.fig)
Saving figure figs/NUC_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/NUC_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/AIT_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/AIT_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/ACC_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/ACC_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/OC_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/OC_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/MOS_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/MOS_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/BC_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/BC_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/MBS_mu_ocean_vs_land_comparison.png...
done.
Saving figure figs/MBS_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/DST01_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/DST02_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/DST03_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/DST04_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/SSLT01_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/SSLT02_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/SSLT03_n_ocean_vs_land_comparison.png...
done.
Saving figure figs/SSLT04_n_ocean_vs_land_comparison.png...
done.
/Users/daniel/workspace/Research/CESM_pcm_activation/aerosol_dists/plot_funcs.py:116: SettingWithCopyWarning: 
A value is trying to be set on a copy of a slice from a DataFrame

See the the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
  df['region'][~np.isnan(df['ocean'])] = 'maritime'
```
