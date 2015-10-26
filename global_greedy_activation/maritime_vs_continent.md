```python
>>> from utils import *
...
>>> %matplotlib inline
>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> sns.set(style='ticks', context='talk')
...
>>> from marc_analysis.utilities import mode_colors
>>> from marc_analysis.analysis import mask_ocean_points
```

Load in the netCDF dataset, and add a mask for data from ocean gridpoints.

```python
>>> data = load_experiment("arg_comp_ensemble/", 'nc')
>>> df = load_experiment("arg_comp_ensemble/", 'csv')
>>> ocean_mask = mask_ocean_points(data)
...
>>> data['region'] = (('index'),
...                   np.where(ocean_mask, 'maritime', 'continental'))
>>> df['region'] = np.where(ocean_mask, 'maritime', 'continental')
Found 10 files
Reading...
Found 10 files
Reading...
    arg_comp_ensemble/single_timestep_ti_0000.csv
    arg_comp_ensemble/single_timestep_ti_0001.csv
    arg_comp_ensemble/single_timestep_ti_0002.csv
    arg_comp_ensemble/single_timestep_ti_0003.csv
    arg_comp_ensemble/single_timestep_ti_0004.csv
    arg_comp_ensemble/single_timestep_ti_0005.csv
    arg_comp_ensemble/single_timestep_ti_0006.csv
    arg_comp_ensemble/single_timestep_ti_0007.csv
    arg_comp_ensemble/single_timestep_ti_0008.csv
    arg_comp_ensemble/single_timestep_ti_0009.csv
```

As a test, how does the distribution of the subgrid vertical velocity change over each region? Are there differences in the vertical velocities that lead to different modes dominating?

```python
>>> ACC = df['MODE_01'] == 'ACC'
>>> df['dominant mode'] = 'other'
>>> df.loc[ACC, 'dominant mode'] = 'ACC'
...
>>> # Remove situations where WSUB is set to the lower limit of 0.2 m/s
... df2 = df.copy()
...
>>> # df2['log(V)'] = np.log10(df2['WSUB'])
... df2.loc[np.isclose(df2['WSUB'], 0.2), 'WSUB'] = np.nan
>>> df2.dropna(inplace=True)
...
>>> mode_colors.update(dict(other='black'))
>>> g = sns.FacetGrid(df2, col="region",
...                   hue="dominant mode", palette=mode_colors,
...                   size=4, aspect=16./12.)
>>> g.map(sns.kdeplot, 'WSUB', shade=True)
>>> g.add_legend()
>>> g.set_ylabels("Probability Density")
>>> # g.set(xlim=(-1.5, 1.0))
... plt.savefig("figs/WSUB_pdf_ocn_vs_lnd.png", dpi=200,
...             bbox_inches='tight')
```

The quick answer would be 'no'; potentially, in maritime regions, there's a preference for slightly weaker updraft speads, and only in the stronger updraft speeds tend to favor other modes dominating. But this is probably covariance with the presence of large sea-salt particles, which we can quickly check -

```python
>>> df2_maritime = df2[df2.region == 'maritime']
>>> g = sns.FacetGrid(df2_maritime,
...                   hue='MODE_01', palette=mode_colors,
...                   size=4, aspect=4./3.)
>>> g.map(sns.kdeplot, 'WSUB', shade=True)
>>> g.add_legend()
>>> g.set_ylabels("Probability Density")
>>> plt.savefig("figs/WSUB_pdf_all_modes_ocn_vs_lnd.png", dpi=200,
...             bbox_inches='tight')
```

![](figs/WSUB_pdf_all_modes_ocn_vs_lnd.png)

In the stronger updrafts, there's a better chance that one of the mixed modes could dominate, but here, sea salt is never the dominant mode.

Going a bit further, does the difference between the activation characteristics of the dominant mode and the full mixture change significantly depending on which mode dominates? Are there differences between the oceanic and continental regimes?

```python
>>> kwargs = dict(shade=True)
...
>>> def _drop_invalid(df):
...     df[~np.isfinite(df)] = np.nan
...     return df.dropna()
...
>>> for region, g in df.groupby('region'):
...
...     # Smax - raw
...     g['dSMAX_01'] = g['SMAX_01'] - g['SMAX_all']
...     df_ACC, df_other = _drop_invalid(np.log10(g.loc[ACC, 'dSMAX_01'])),\
...                        _drop_invalid(np.log10(g.loc[~ACC, 'dSMAX_01']))
...
...     fig = plt.figure(figsize=(15, 3.5))
...     ax = fig.add_subplot(131)
...     kwargs['ax'] = ax
...
...     sns.kdeplot(df_ACC, label="ACC dominant", color=mode_colors['ACC'],
...                 **kwargs)
...     sns.kdeplot(df_other, label="Other dominant", color='k', **kwargs)
...     sns.despine(fig)
...     ax.set_xlabel("$\Delta$ log$_{10}$(Smax)")
...     ax.set_ylabel("density function")
...     # ax.set_xlim(0, 0.01)
...
...     # Smax - percent
...     df_ACC = _drop_invalid(100.*g.loc[ACC, 'dSMAX_01']/
...                             g.loc[ACC, 'SMAX_all'])
...     df_other = _drop_invalid(100.*g.loc[~ACC, 'dSMAX_01']/
...                              g.loc[~ACC, 'SMAX_all'])
...
...     ax = fig.add_subplot(132)
...     kwargs['ax'] = ax
...     kwargs['legend'] = False
...
...     sns.kdeplot(df_ACC, label="ACC dominant", color=mode_colors['ACC'],
...                 **kwargs)
...     sns.kdeplot(df_other, label="Other dominant", color='k', **kwargs)
...     sns.despine(fig)
...     ax.set_xlabel("$\Delta$ Smax (%)")
...     ax.set_ylabel("density function")
...     ax.set_xlim(0)
...
...     # As a reference, also plot the distribution of Smax (1)
...     df_ACC, df_other = _drop_invalid(np.log10(g.loc[ACC, 'SMAX_01'])), \
...                        _drop_invalid(np.log10(g.loc[~ACC, 'SMAX_01']))
...
...     ax = fig.add_subplot(133)
...     kwargs['ax'] = ax
...
...     sns.kdeplot(df_ACC, label="ACC dominant", color=mode_colors['ACC'],
...                 **kwargs)
...     sns.kdeplot(df_other, label="Other dominant", color='k', **kwargs)
...     sns.despine(fig)
...     ax.set_xlabel("log$_{10}$(Smax) (1st mode)")
...     ax.set_ylabel("density function")
...
>>> #     fig.suptitle(region, fontweight='bold', fontsize=16)
...     plt.savefig("figs/delta_Smax_lnd_vs_ocn_{}.png".format(region),
...                 dpi=200, bbox_inches='tight')
/Users/daniel/anaconda/envs/marc_analysis/lib/python3.4/site-packages/ipykernel/__main__.py:10: SettingWithCopyWarning: 
A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead

See the the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
```

![](figs/delta_Smax_lnd_vs_ocn_continental.png)
![](figs/delta_Smax_lnd_vs_ocn_maritime.png)

There is very clearly a difference between the characteristics depending on which mode domiantes. When ACC dominates, $\Delta$ Smax tends to be much smaller. Put in other terms, when ACC dominates, the activation characteristics of the full aerosol mixture tend to closely follow that when only ACC is present. This is true in either maritime or continental regimes.

To illustrate this point in more detail, we can compare the distributions of Smax for each region:

```python
>>> df['log(Smax)'] = np.log10(df['SMAX_01'])
...
>>> g = sns.FacetGrid(df, hue='region',
...                   palette=dict(continental='b', maritime='g'),
...                   size=4, aspect=4./3., legend_out=False)
>>> g.map(sns.kdeplot, 'log(Smax)', shade=True)
>>> g.add_legend()
>>> g.set_ylabels("Probability Density")
>>> plt.savefig("figs/smax_lnd_v_ocn.png", dpi=200, bbox_inches='tight')
```

![Smax comparison - land vs ocean](figs/smax_lnd_v_ocn.png)

Not much of a difference; if anything, in the continental regimes, there is a slightly wider distribution (more variability). But they're not grossly different.
