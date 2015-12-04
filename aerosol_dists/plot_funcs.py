"""
Plotting utility functions for the comparison box plots.
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from collections import OrderedDict

mode_colors = {
    # Sulfate - red
    'NUC': '#fee0d2',
    'AIT': '#fc9272',
    'ACC': '#ef3b2c',
    # Organic - green
    'OC': '#4daf4a',
    'MOS': '#ffff33', # green + red = yellow
    # Black - blue,
    'BC': '#377eb8',
    'MBS': '#984ea3', # blue + red = purple
    # Dust - shades of brown
    'DST01': '#f6e8c3', 'DST02': '#dfc27d',
    'DST03': '#bf812d', 'DST04': '#8c510a',
    # Sea Salt - shades of teal
    'SSLT01': '#c7eae5', 'SSLT02': '#80cdc1',
    'SSLT03': '#35978f', 'SSLT04': '#01665e',
}
all_modes = ['NUC', 'AIT', 'ACC', 'OC', 'MOS', 'BC', 'MBS',
             'DST01', 'DST02', 'DST03', 'DST04',
             'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04', ]

# Sub-sampling data logic
levels = {
    'strat': slice(30, 300),
    'near-sfc': slice(850, 1100),
    'trop': slice(500, 850),
}

lats = {
    'tropics': slice(-20, 20),
    'nh': slice(20, 80),
    'sh': slice(-80, -20),
    'global': slice(-80, 80),
}

time_subset = slice("2001-10", None)

def compare_boxplots(mode, param, lat_lev_generator, data):

    sns.set(style='ticks', context="talk")

    level_labels = set()
    lat_labels = set()
    comp_data = []

    for subset in lat_lev_generator():
        (lev, lev_slice), (lat, lat_slice) = subset

        level_labels.add(lev)
        lat_labels.add(lat)

        data_subset = data.sel(lev=lev_slice, lat=lat_slice)
        data_subset = data_subset.isel(time=-1, nbnd=0)

        ## Original promotion to DataFrame
        # n_aer = data_subset[param + mode]
        # df = pd.DataFrame({mode: n_aer.values.ravel(),
        #                    'lev': lev,
        #                    'lat': lat, })
        ## Automatic promotion to DataFrame using xray machinery
        df = data_subset[[param + mode, ]].to_dataframe()
        df.rename(columns={param+mode: mode}, inplace=True)
        df.reset_index(inplace=True) # MultiIndex -> Normal index with
                                     # dimension values as columns
        df['lat'] = lat
        df['lev'] = lev

        comp_data.append(df)
    all_df = pd.concat(comp_data, ignore_index=True)

    g = sns.factorplot(x="lat", y=mode, hue="lev", data=all_df,
                       hue_order=['near-sfc', 'trop', 'strat'],
                       order=['global', 'nh', 'tropics', 'sh'],
                       kind='box', aspect=16./10., size=5,
                       whis=[10., 90.], color=mode_colors[mode],
                       showfliers=False)

    plt.semilogy()
    if param == 'mu':
        plt.ylim(1e-4, 1.)
        plt.ylabel("Mean Size, {} [micron]".format(mode))
    elif param == 'n':
        plt.ylim(1e-2, 1e4)
        plt.ylabel("Number Conc, {} [cm^-3]".format(mode))

    return g

def compare_maritime_vs_continent(mode, param, data):

    comp_data = []
    lat_labels = [lat for lat, _ in lats.items()]
    lev_slice = levels['near-sfc']
    for lat, lat_slice in lats.items():
        data_subset = (data
                           # .isel(time=(slice(-10, -1)), nbnd=0)
                           .isel(time=-1, nbnd=0)
                           .sel(lat=lat_slice, lev=lev_slice))

        # Basd on the ocean mask, add a categorical field indicating
        # the region
        df = data_subset[[param+mode, 'ocean']].to_dataframe()
        df.reset_index(inplace=True)
        df.rename(columns={param+mode: mode}, inplace=True)
        df['region'] = 'continent'
        # NaN -> continent cell; masking with ~np.isnan(ocean) chooses
        # the points which aren't NaN, and thus are ocean cells
        df['region'][~np.isnan(df['ocean'])] = 'maritime'
        del df['ocean']
        df['lat'] = lat

        comp_data.append(df)
    all_df = pd.concat(comp_data, ignore_index=True)

    g = sns.factorplot(x='lat', y=mode, hue='region', data=all_df,
                       palette=dict(maritime='b', continent='grey'),
                       kind='box', aspect=16./10., size=5, whis=[10., 90.],
                       showfliers=False)

    plt.semilogy()
    if param == 'mu':
        plt.ylim(1e-4, 1.)
        plt.ylabel("Mean Size, {} [micron]".format(mode))
    elif param == 'n':
        plt.ylim(1e-2, 1e4)
        plt.ylabel("Number Conc, {} [cm^-3]".format(mode))

    return g

def _to_str(x):
    if int(x) == x:
        return str(int(x))
    else:
        p = int(np.floor(np.abs(np.log10(x))))
        fmt_str = "{val:.{p:d}f}".format(val=x, p=p)
        return fmt_str
def mode_percentiles(data, param='n', percentiles=[.1, 1., 2., 5., 10.]):
    """ Compute percentiles of a specific parameter for each mode in a given
    dataset. 'percentiles' will be mirrored around the median, so only the
    lower values need to be provided. """

    all_mode_percentiles = {}

    for mode in all_modes:

        try:
            aer_data = data[param + mode].values.ravel()
        except KeyError:
            continue

        pcts = np.array(percentiles + [100. - p for p in percentiles[::-1]])
        tiles = list(np.percentile(aer_data, pcts))

        vals = zip(['min', ] + [_to_str(p) for p in pcts] + ['max', ],
                   [np.min(aer_data), ] + tiles + [np.max(aer_data), ])

        mode_percentiles = pd.Series(OrderedDict(vals))
        all_mode_percentiles[mode] = mode_percentiles

    return pd.DataFrame(all_mode_percentiles)


def aerosol_boxplots(mode, data, percentile=10.,
                     size=5., aspect=16./10.,
                     ratio=5):
    """ Generated paired boxplots and distribution plots of the
    number and size parameters for a given aerosol distribution. """

    sns.set(style='ticks', context="talk")

    n_aer = data['n' + mode].values.ravel()
    l_pct_N = np.percentile(n_aer, percentile)
    r_pct_N = np.percentile(n_aer, 100.-percentile)

    skip_mu = True
    try:
        mu_aer = data['mu' + mode].values.ravel()
        l_pct_mu = np.percentile(mu_aer, percentile)
        r_pct_mu = np.percentile(mu_aer, 100.-percentile)
        skip_mu = False
    except KeyError:
        mu_aer = [1., ]
        l_pct_mu = r_pct_mu = 1.
    color = mode_colors[mode]

    # Set up the figure grid:
    #
    #    KDE - Number  |   KDE - Size
    # -----------------|-----------------
    #                  |
    #    Box - Number  |   Box - Size
    #                  |
    fig = plt.figure(figsize=(2*size, size))
    gs = plt.GridSpec(ratio+1, 2*ratio)

    box_N_ax = fig.add_subplot(gs[1:, :ratio])
    dist_N_ax = fig.add_subplot(gs[0, :ratio], sharex=box_N_ax)
    box_mu_ax = fig.add_subplot(gs[1:, ratio:])
    dist_mu_ax = fig.add_subplot(gs[0, ratio:], sharex=box_mu_ax)

    sns.boxplot(np.log10(n_aer), whis=[10., 90.], showfliers=False,
                ax=box_N_ax, color=color)
    sns.kdeplot(np.log10(n_aer), clip=np.log10([l_pct_N, r_pct_N]),
                ax=dist_N_ax, color=color, shade=True)

    if not skip_mu:
        sns.boxplot(np.log10(mu_aer), whis=[10., 90.], showfliers=False,
                    ax=box_mu_ax, color=color)
        sns.kdeplot(np.log10(mu_aer), clip=np.log10([l_pct_mu, r_pct_mu]),
                    ax=dist_mu_ax, color=color, shade=True)

    # Set consistent plot limits
    box_N_ax.set_xlim(-3, 4)
    box_mu_ax.set_xlim(-3, 0)

    # Turn off tick visibility on x-axis of distribution plots
    plt.setp(dist_N_ax.get_xticklabels(), visible=False)
    plt.setp(dist_mu_ax.get_xticklabels(), visible=False)

    plt.setp(dist_N_ax.get_yticklabels(), visible=False)
    plt.setp(dist_mu_ax.get_yticklabels(), visible=False)
    plt.setp(box_N_ax.get_yticklabels(), visible=False)

    # Turn off ticks on density axis for boxplots
    plt.setp(dist_N_ax.yaxis.get_majorticklines(), visible=False)
    plt.setp(dist_mu_ax.yaxis.get_majorticklines(), visible=False)
    plt.setp(dist_N_ax.xaxis.get_majorticklines(), visible=False)
    plt.setp(dist_mu_ax.xaxis.get_majorticklines(), visible=False)

    # Axis labels
    box_N_ax.set_ylabel(mode)
    box_N_ax.set_xlabel("log10(Number Concentration [cm-3])")
    box_mu_ax.set_xlabel("log10(Mean Radius [micron])")

    sns.despine(fig=fig, trim=True)
    sns.despine(ax=dist_N_ax, bottom=True, left=True)
    sns.despine(ax=dist_mu_ax, bottom=True, left=True)
    sns.despine(ax=box_mu_ax, left=True)
    fig.tight_layout()
    fig.subplots_adjust(hspace=.2, wspace=1.)

    return fig

