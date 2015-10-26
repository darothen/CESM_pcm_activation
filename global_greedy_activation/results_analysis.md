```{.python .input  n=1}
import glob, os

%matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks')

import numpy as np
import pandas as pd
import xray

from utils import *
```

Load in some of the partial data from the actual simulation. We'll load the CSV
data for simplicity, and find them all by using the `glob` module.

```{.python .input  n=2}
DATA_DIR = "arg_comp_ensemble/"
df = load_experiment(DATA_DIR, 'csv')
ds = load_experiment(DATA_DIR, 'nc')
```

```{.json .output n=2}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "Found 10 files\nReading...\n    arg_comp_ensemble/single_timestep_ti_0000.csv\n    arg_comp_ensemble/single_timestep_ti_0001.csv\n    arg_comp_ensemble/single_timestep_ti_0002.csv\n    arg_comp_ensemble/single_timestep_ti_0003.csv\n    arg_comp_ensemble/single_timestep_ti_0004.csv\n    arg_comp_ensemble/single_timestep_ti_0005.csv\n    arg_comp_ensemble/single_timestep_ti_0006.csv\n    arg_comp_ensemble/single_timestep_ti_0007.csv\n    arg_comp_ensemble/single_timestep_ti_0008.csv\n    arg_comp_ensemble/single_timestep_ti_0009.csv\nFound 10 files\nReading...\n"
 }
]
```

What proportion of the modes correspond to the dominant activating mode at each
iteration? We can do some summary statistics on this:

**First mode**

```{.python .input  n=3}
grouped_i1 = df.groupby('MODE_01')
n_i1 = grouped_i1.size()
pct_i1 = 100.*n_i1/len(df)
df_i1 = pd.DataFrame({'total': n_i1, 'pct': pct_i1})

print (df_i1)
```

```{.json .output n=3}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "               pct  total\nMODE_01                  \nACC      96.495970  57473\nDST01     0.142713     85\nMBS       1.145064    682\nMOS       2.216253   1320\n"
 }
]
```

**Second mode**

```{.python .input  n=4}
grouped_i2 = df.groupby(['MODE_01', 'MODE_02'])
n_i2 = grouped_i2.size() # this is a DataFrame!

# group by the first iteration
pct_i2 = n_i2.groupby(level=0) \
             .apply(lambda g: 100.*g/g.sum()) # compute the fraction of
                                              # of 2nd iter mode contribs
                                              # given that first iter mode

df_i2 = pd.DataFrame({'total': n_i2, 'pct': pct_i2})

print (df_i2)
```

```{.json .output n=4}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "                       pct  total\nMODE_01 MODE_02                  \nACC     DST01     0.742958    427\n        MBS      10.321368   5932\n        MOS      36.164808  20785\n        SSLT01   52.770866  30329\nDST01   ACC      32.941176     28\n        DST02    44.705882     38\n        MOS      22.352941     19\nMBS     ACC      36.070381    246\n        DST01     0.586510      4\n        MOS      63.343109    432\nMOS     ACC      45.681818    603\n        DST01     1.590909     21\n        MBS      52.727273    696\n"
 }
]
```

### So very clearly, the **ACC** mode is really important; 95% of the time it's
the dominant mode. How often is it either the 1st- or 2nd-most important mode?

```{.python .input}
ACC_tot = n_i1.ix['ACC'] + n_i2.xs('ACC', level='MODE_02').sum()
#       ACC iter 1 total      n from iter 2 where MODE_02 == ACC
# alternatively: use `sel` from below
print( 100.*ACC_tot/n_i1.sum())
```

... 99% of the time!

```{.python .input}
# Export mode data as JSON for sunburst viz

# 1) Construct the hierarchical dataset
from itertools import permutations
mode_names = ['AIT', 'ACC', 'MBS', 'MOS', 'DST01', 'DST02', 'SSLT01']
multi_labels = list(permutations(mode_names, 3))

full_index = pd.MultiIndex.from_tuples(multi_labels, 
                                       names=['MODE_01', 
                                              'MODE_02', 
                                              'MODE_03'])

data = df.groupby(['MODE_01', 'MODE_02', 'MODE_03']).size()
# PCT by LEVEL
# data_pct = data.groupby(level=0) \
#                .apply(lambda g: 100.*g/g.sum()) # compute the fraction of
#                                                 # of 2nd iter mode contribs
#                                                 # given that first iter mode
# PCT by ALL
data_pct = 100.*data/(np.sum(data))

data = pd.DataFrame({'pct': data_pct, 'count': data})
data_full = data.reindex(full_index, fill_value=0, copy=True)
mask = data_full['pct'] > 0.0

data_out = data_full[mask]
print (data_out)
```

```{.python .input}
from collections import OrderedDict as od

# Create JSON to write to disk
m1_children = []
m1s = set(data_out.index.get_level_values(0))
for m1 in sorted(m1s):
    m2s = set(data_out.ix[m1].index.get_level_values(0))

    m2_children = []
    for m2 in sorted(m2s):
        m3s = set(data_out.ix[m1,m2].index.get_level_values(0))                
        m3_children = [ { 'name': m3,
                          'size': int(data_out.ix[m1,m2,m3]['count']) } 
                       for m3 in sorted(m3s) ]
        
        m2_children.append( { 'name': m2, 'children': m3_children })
    m1_children.append( { 'name': m1, 'children': m2_children } )

data_json = {'name': 'mode_data', 'children': m1_children}

import json
with open("mode.json", "w") as f:
    json.dump(data_json, f)
    
import pprint
pprint.pprint(data_json)
```

Is there something special about the situations where **ACC** *isn't* one of the
two most important modes? We can split the original dataset into two segments
and visualize using a `factorplot` on that identifier.

```{.python .input}
#sel = (df['MODE_01'] == 'ACC') | (df['MODE_02'] == 'ACC')
sel = df['MODE_01'] == 'ACC'
yes_ACC = df[sel]
no_ACC = df[~sel]

# Label the samples based on of ACC was in the first two iters
df['ACC_important'] = True
df.loc[~sel, 'ACC_important'] = False

# Plot KDEs of each of the parameters, split along ACC_important
keys = [k for k in df.keys() if not "_" in k]
n_plots = len(keys)
n_cols = 4
n_rows = int(n_plots / n_cols)

fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*3.5, n_rows*2))
for i, (key, ax) in enumerate(zip(keys, axs.ravel())):
    
    kde_args = dict(ax=ax, shade=True)
    if i > 0: kde_args['legend'] = False
    
    df_ACC, df_other = df.loc[sel, key], df.loc[~sel, key]
    if key.startswith("n") or key.startswith("mu"):
        df_ACC = np.log10(df_ACC)
        df_other = np.log10(df_other)
        key = "log10(%s)" % key
    try:
        sns.kdeplot(df_ACC, label="ACC dominant", 
                    color='k', **kde_args)
        sns.kdeplot(df_other, label="Other dominant",
                    color='r', **kde_args)
    except:
        print('something broke')
    
    ax.set_xlabel(key)
    plt.setp(ax.get_xticklabels(), rotation=20)
    #labels = [t.get_text() for t in ax.get_xticklabels()]
    #ax.set_xticklabels(labels, rotation=45)
    
    if (i % n_cols) == 0:
        ax.set_ylabel("density function")

sns.despine(fig=fig)
plt.tight_layout()
```

A couple things stand out, ranked form *curious* to *obvious* -

1. **ACC** dominates more often in the SH than the NH, but the *opposite* is
true for when other modes dominate. Does this have something to do with inter-
hemispheric differences in the aerosol size distributions and other factors
identified below?

2. Generally speaking, when **other** modes dominate, ACC there are (a) *fewer*
ACC particles, and (b) *smaller* ACC particles. This is likely because the size
distribution moments have a strong positive correlation?

3. The opposite is true for MBS/MOS. When ACC doesn't dominate, they tend to be
*smaller* but *more numerous*. Do the same correlations between moments exit for
these two modes?

The next few analyses dig into these questions.

## 1) Splitting by hemisphere

```{.python .input}
sel_sh = df['lat'] < 0.
df['hemisphere'] = 'NH'
df.loc[sel_sh, 'hemisphere'] = 'SH'

n_hemi = df.groupby('hemisphere').size()
print (n_hemi)

for key in keys:
    if key in ['lat', 'lev', 'lon', 'time']: continue
    
    df_trunc = df[['hemisphere', 'ACC_important', key]]
    if key.startswith("n") or key.startswith("mu"):
        df_trunc[key] = np.log10(df_trunc[key])
        label = "log10(%s)" % key
    else: label = key
    
    print(key)
    with sns.color_palette(['k', 'r']):
        g = sns.FacetGrid(df, col='hemisphere', hue='ACC_important', legend_out=False,
                          size=2.5, aspect=3.5/2.5, hue_order=[True, False])
        g.map(sns.kdeplot, key, shade=True)
        g.add_legend(title="ACC dominant?")
        g.set_ylabels("density function")
        g.set_xlabels(label)
        g.set_xticklabels(rotation=20)
        
```

Splitting by hemisphere doesn't really change anything; there aren't any wildly
different distributions in the input parameters and the ensuing dominant mode.
Obviously in the NH, we prove colder temperatures. What's curious is the
preponderance of SH vs NH samples; why should there be more valid samples from
the SH than the NH? This is worth considering a bit more. First, let's try to
plot some sort of heatmap of where the parameter samples came from.

```{.python .input}
# Jitter the lat/lon coordinates to represent sampling from gridboxes,
# which in these simulations were 2.5 x 1.9 degrees
dlon, dlat = 2.5, 1.9
npoints = len(df)
df_map = df.copy()

def rand_fact(npoints, width):
    """ Rescale a number from [0, 1) to [-width, width) """
    points = 2.*np.random.random(npoints) - 1. 
    points *= width
    return points

l1 = df_map['lon'].copy()
l2 = df_map['lat'].copy()

df_map['lon'] = df_map['lon'] + rand_fact(npoints, dlon/2.)
df_map['lat'] = df_map['lat'] + rand_fact(npoints, dlat/2.)

lon, lat = df_map.lon, df_map.lat
# # Correct lon: 
# # 1) some values may be < 0 or > 360, map these into [0, 360]
# lon[lon < 0] = lon[lon < 0] + 360.
# lon[lon > 360] = lon[lon > 360] - 360.
# # 2) map from [0, 360] -> [-180, 180]
# lon -= 180.

df_map['lon'] = lon[:]
```

```{.python .input}
import cartopy.crs as ccrs
import cartopy.feature

lon, lat = df_map.lon, df_map.lat
# Correct lon: 
# 1) some values may be < 0 or > 360, map these into [0, 360]
# lon[lon < 0] = lon[lon < 0] + 360.
# lon[lon > 360] = lon[lon > 360] - 360.
# # 2) map from [0, 360] -> [-180, 180]
# lon -= 180.

proj = ccrs.PlateCarree()
cmap = sns.light_palette("navy", 12, as_cmap=True)

with sns.color_palette("Reds"):
    ax = plt.axes(projection=proj)
    hb = ax.hexbin(lon, lat, gridsize=(50, 15), bins='log', 
                   transform=proj, cmap=cmap)
    ax.coastlines()
    ax.set_global()
```

It's not immediately obvious why we'd have see this separation between the two
hemispheres. However, given that we're over-sampling remote maritime conditions,
an obvious hypothesis would be bias in the $\mu$-$N$ relationships for the
sulfate and mixed modes. If increasing $N$ tends to decrease $\mu$ below our
cutoff threshholds, then we could be in for some trouble! We should explicitly
plot those relationships for the four modes, segregating by hemisphere.

```{.python .input}
## Regressions
modes = ['AIT', 'ACC', "MOS", "MBS"]
df_aero = df[['mu%s' % m for m in modes] + ['n%s' % m for m in modes]] 
df_aero = np.log10(df_aero)
df_aero['hemisphere'] = df['hemisphere']
df_aero['ACC_important'] = df['ACC_important']

for mode in modes:
    with sns.color_palette(['k', 'r']):
        g = sns.lmplot('mu'+mode, 'n'+mode, df_aero, 
                       row='ACC_important', col='hemisphere', hue='ACC_important',
                       size=3., aspect=1.5, 
                       markers='.', scatter_kws=dict(alpha=0.7))
        
        for ax in g.axes.ravel():
            ax.set_ylim(0, 4)
            ax.set_xlim(-3, 0)
```

```{.python .input}
## Distributions
cmap = sns.light_palette("navy", 12, as_cmap=True)

for mode in modes: 
    g = sns.FacetGrid(df_aero, row='ACC_important', col="hemisphere",
                      size=3, aspect=1.5, margin_titles=True)
    g.map(plt.hexbin, 'n'+mode, 'mu'+mode, gridsize=(24, 12), edgecolor='none',
          bins='log', cmap=cmap)
```

So there are definitely some non-trivial differences between the aerosol in the
NH/SH. But it's probably better to focus on the whole rather than break it into
the parts that I've done here.

---

## $\Delta$Smax, $\Delta$Nact

Is $\Delta$Smax impacted by whether or not ACC is important?

```{.python .input}
fig, axs = plt.subplots(1, 4, figsize=(15, 3))
for i in xrange(1, 5):
    i_str = "%02d" % i
    df['dSMAX_' + i_str] = df['SMAX_' + i_str] - df['SMAX_all']
    ax = axs.ravel()[i-1]
    sns.kdeplot(df['dSMAX_%02d' % i], ax=ax, shade=True, color='k',
                label=None)
sns.despine(fig)
```

```{.python .input}
sel = df['MODE_01'] == 'ACC'

kwargs = dict(shade=True)

# Smax - raw
df_ACC, df_other = df.loc[sel, 'dSMAX_01'], df.loc[~sel, 'dSMAX_01']

fig = plt.figure(figsize=(15, 3.5))
ax = fig.add_subplot(131)
kwargs['ax'] = ax

sns.kdeplot(df_ACC, label="ACC dominant", color='k', **kwargs)
sns.kdeplot(df_other, label="Other dominant", color='r', **kwargs)
sns.despine(fig)
ax.set_xlabel("$\Delta$ Smax")
ax.set_ylabel("density function")
ax.set_xlim(0)

# Smax - percent
df_ACC = 100.*df.loc[sel, 'dSMAX_01']/df.loc[sel, 'SMAX_all']
df_other = 100.*df.loc[~sel, 'dSMAX_01']/df.loc[~sel, 'SMAX_all']

ax = fig.add_subplot(132)
kwargs['ax'] = ax
kwargs['legend'] = False

sns.kdeplot(df_ACC, label="ACC dominant", color='k', **kwargs)
sns.kdeplot(df_other, label="Other dominant", color='r', **kwargs)
sns.despine(fig)
ax.set_xlabel("$\Delta$ Smax (%)")
ax.set_ylabel("density function")
ax.set_xlim(0)

# As a reference, also plot the distribution of Smax (1)
df_ACC, df_other = np.log10(df.loc[sel, 'SMAX_01']), np.log10(df.loc[~sel, 'SMAX_01'])

ax = fig.add_subplot(133)
kwargs['ax'] = ax

sns.kdeplot(df_ACC, label="ACC dominant", color='k', **kwargs)
sns.kdeplot(df_other, label="Other dominant", color='r', **kwargs)
sns.despine(fig)
ax.set_xlabel("log$_{10}$(Smax) (1st mode)")
ax.set_ylabel("density function")

```

```{.python .input}
from bokeh import mpl
from bokeh.plotting import show, output_notebook
from IPython import display
```

```{.python .input}
r = df.ix[0]
keys = ['SMAX_01', 'SMAX_02', 'SMAX_03', 'SMAX_04']
idx = np.array(range(1, len(keys)+1))

r_smax = r[keys]*100.
plt.plot(idx, r_smax, marker='s', color='r')
plt.hlines(r['SMAX_all']*100., idx[0], idx[-1],
           alpha=0.5, linestyle='--')
plt.xticks(idx)
plt.xlim(idx[0]-0.1, idx[-1]+0.1)

plt.ylabel("Smax (%)")
plt.xlabel("iteration")
sns.despine(offset=20)

```
