#!/usr/bin/env python
""" Plot one-one comparisons between parcel model and
a given activation scheme for a specific sampling
experiment. """

import os, pickle
import numpy as np
import pandas as pd
import sklearn.metrics as skm
from statsmodels.distributions import ECDF

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set(style='ticks', context='talk')

def get_data(exp_name, data_root="data/", tidy=False):
    """ Load a DataFrame with the raw sampling results from a given
    experiment. """
    path_to_data = os.path.join(data_root, exp_name)
    which = "tidy" if tidy else "results"
    results_data_fn = "{}_sampling_{}.csv".format(exp_name, which)

    # Read in the sampling results
    results_df = pd.read_csv(
        os.path.join(path_to_data, results_data_fn), index_col=0
    )

    return results_df


def get_stats(exp_name, result='Smax',
              scaling='normal', override=False):
    """ Retrieve the pre-computed statistics from the sampling experiments.

    Parameters
    ----------
    result : str
        'Smax', 'Neq', 'Nderiv_Neq', etc.
    scaling : str
        'normal' or 'log10'

    """

    stats_fn = "{}_stats_processed.p".format(exp_name)
    stats = pd.read_pickle(stats_fn)

    if override:
        return stats

    s = stats.copy()
    # Choose variable
    s = s.T.xs(result, level='result').T
    # Choose sampling
    s = s.T.xs(scaling, level='scaling').T

    return s


def summary_stats(obs, act):
    """ Create a Series with summary statistics comparing two
    array-like datasets. """

    mae = skm.mean_absolute_error(act, obs)
    r2 = skm.r2_score(act, obs)
    rmse = np.sqrt(np.sum((obs-act)**2.)/len(act))
    nrmse = rmse/np.sqrt(np.sum((act**2.)/len(act)))

    rel_err = 100.*(obs - act)/act
    # Mask egregiously high values (1000% error) which screw up the spread
    rel_err = rel_err[np.abs(rel_err) <= 1000.]
    mre = np.mean(rel_err)
    mre_std = np.std(rel_err)

    stats = pd.Series({
        'mae': mae, 'r2': r2, 'rmse': rmse, 'nrmse': nrmse,
        'mre': mre, 'mre_std': mre_std,
    })

    return stats


def clean_df(df, lower, upper):
    """ Clean a DataFrame by dropping values outside the range
    [lower, upper]

    Parameters
    ----------
    df : DataFrame
        The DataFrame to clip
    lower, upper : floats
        Lower and upper bound of valid data range

    """

    df[df > upper] = np.nan
    df[df < lower] = np.nan
    df = (
        df
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )

    return df


def plot_oneone(parcel_data, param_data, var_label, param_label,
                color_data=None, color_lognorm=True,
                color_ticks=None, color_bar=False, ax=None,
                lims=(1e-3, 10.), loglog=True, scatter_kws=None):
    """ Create a one-one plot comparing two different data sources.

    Parameters
    ----------
    parcel_data, param_data : array-like
        The x- and y-axis datasets, respectively
    var_label : str
        Label to be printed on x- and y-axes indicating variable name and
        information
    param_label : str
        Label to be printed on y-axis indicating parameterization name
    color_data : array-like
        Data to use for coloring points on the plot
    color_lognorm : logical
        Lognormalize the colorbar?
    color_ticks : array-like
        Ticks to label on color bar, if being included by `color_bar`
    color_bar : logical
        Include color bar indicating coloring?
    lims : length-2 tuple
        x/y-axis limits
    loglog : logical (default = True)
        Automatically log-log both axes
    scatter_kws : dict
        Additional keyword arguments to pass to `scatter` function.

    Returns
    -------
    fig, ax
        The figure and axis for drawing, respectively.

    """

    if ax is None:
        fig = plt.figure(figsize=(5., 5./(16./10.)))
        ax = fig.add_subplot(111)
    else:
        plt.sca(ax)
        fig = plt.gcf()

    base_scatter_kws = dict(marker='.', s=30, edgecolor='none', alpha=0.8,
                            cmap=plt.get_cmap("viridis"))
    # Update plot settings with user-defined parameters
    if scatter_kws is not None:
        base_scatter_kws.update(scatter_kws)
    # Add coloring if available
    if color_data is not None:
        base_scatter_kws['c'] = color_data
        if color_lognorm:
            base_scatter_kws['norm'] = mpl.colors.LogNorm()
        if color_ticks is not None:
            base_scatter_kws['vmax'] = color_ticks[-1]

    scatr = ax.scatter(parcel_data, param_data,
                       **base_scatter_kws)

    oo = np.linspace(*lims, num=100)
    ax.plot(oo, oo, color='grey', lw=3)
    ax.plot(oo, oo*0.5, color='k', lw=1, alpha=0.8)
    ax.plot(oo, oo*1.5, color='k', lw=1, alpha=0.8)

    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    if loglog:
        ax.loglog()

    # Set tick formatter to 1 significant digit
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%1g"))
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%1g"))

    ax.set_xlabel("{} - parcel model".format(var_label))
    ax.set_ylabel("{} - {}".format(var_label, param_label))

    # Embed a colorbar in the left-hand subplot
    if color_bar:
        cb_bb = mpl.transforms.Bbox.from_bounds(0.75, 0.05, 0.025, 0.4)
        trans = ax.transAxes + fig.transFigure.inverted()
        l, b, w, h = mpl.transforms.TransformedBbox(cb_bb, trans).bounds
        cb_axes = fig.add_axes([l, b, w, h])
        cb = plt.colorbar(scatr, cax=cb_axes,
                          ticks=color_ticks, format="%1g",
                          orientation='vertical', )
        cb.ax.tick_params(labelsize=10)
        cb.outline.set_visible(False)

    sns.despine(ax=ax)

    return fig, ax


def plot_cdf(data, line_kws=None,
             pct_levs=np.linspace(0., 100., 21), ax=None):

    # fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    # plt.subplots_adjust(wspace=0.25, bottom=0.15)
    # ax_cdf, ax_pdf = axs
    if ax is None:
        fig = plt.figure(figsize=(5., 5./(16./10.)))
        ax = fig.add_subplot(111)
    else:
        plt.sca(ax)
        fig = plt.gcf()

    # Compute empirical CDFs
    # data_cdf = ECDF(data.ravel())
    data_percentiles = [np.percentile(data.ravel(), x)
                        for x in pct_levs]

    default_line_kws = dict(color='k', lw=5)
    if line_kws is not None:
        default_line_kws.update(line_kws)
    ax.plot(data_percentiles, pct_levs/100., **default_line_kws)

    ax.set_ylim(0, 1)
    ax.set_ylabel("Cumulative Probability")

    return fig, ax

    # # PDFs
    # plot_kwargs = {}
    # if loglog:
    #     plot_kwargs['hist'] = True
    #     plot_kwargs['kde'] = False
    #     plot_kwargs['bins'] = np.logspace(*[np.log10(l) for l in lims])
    # else:
    #     plot_kwargs['hist'] = False
    #     plot_kwargs['kde'] = True
    #     plot_kwargs['bins'] = np.linspace(*lims)
    # ax_pdf = sns.distplot(parcel_data,
    #                       color='k', label="parcel model", ax=ax_pdf,
    #                       kde_kws={'lw': 5},
    #                       **plot_kwargs)
    # ax_pdf = sns.distplot(data_df['scheme'],
    #                       label=param_name, ax=ax_pdf,
    #                       kde_kws={'linestyle': 'dashed'},
    #                       **plot_kwargs)
    # ax_pdf.set_xlim(*lims)
    # ax_pdf.set_xlabel(var_name)
    # if loglog:
    #     ax_pdf.set_xscale('log')
    # ax_pdf.set_ylabel("Probability Density")
    # ax_pdf.legend(loc="best")
