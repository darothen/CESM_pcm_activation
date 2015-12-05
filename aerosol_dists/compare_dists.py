#!/usr/bin/env python
"""
Compare aerosol distribution parameters from CESM/MARC output.

This script will produce boxplots of the distributions of aerosol size
distribution parameters from a specified dataset. It digs further by
subsetting variables based on maritime/land, different latitude ranges,
and different altitudes.

"""
from __future__ import print_function, division

import os
import numpy as np
import pandas as pd
import xray

from itertools import product
from scipy.stats import percentileofscore

import plot_funcs
from plot_funcs import (lats, levels, all_modes, aerosol_boxplots,
                        all_modes, mode_percentiles)
cb = plot_funcs.compare_boxplots
cmc = plot_funcs.compare_maritime_vs_continent

from marc_analysis import save_figure
from marc_analysis.analysis import extract_feature

import matplotlib.pyplot as plt
plt.ioff()
import seaborn as sns
sns.set(style='ticks', context="talk")

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("aerosol_ds", type=str, metavar="aerosol_data.nc",
                    help="Input file containing the aerosol size"
                         " distribution parameters")
parser.add_argument('-m', "--mode", nargs="+", metavar='',
                    choices=all_modes, default=all_modes,
                    help="Mode(s) to run through plots")
parser.add_argument("-o", "--out_path", type=str,
                    default="/Users/daniel/Dropbox_MIT/Research/results/CESM_pcm_activation",
                    help="Path to save figures")

def all_subsets(supersets=True):
    """ Generate all subsets of a given set of combinations. """
    for lev, lev_slice in levels.items():
        if supersets:
            yield ((lev, lev_slice), ("global", slice(-80, 80)))

        for lat, lat_slice in lats.items():
            if lat == 'global':
                continue
            yield (lev, lev_slice), (lat, lat_slice)

if __name__ == "__main__":

    args = parser.parse_args()

    # Convenience function for prepending output path
    # _out_path = lambda s: os.path.join(args.out_path, s)
    _out_path = lambda s: s

    # Nudge times to the year 2000
    data = xray.open_dataset(args.aerosol_ds,
                             decode_times=False)
    times = data.coords.to_dataset().time
    times += 2000.*365
    data = xray.decode_cf(data)

    # Global troposphere slice for quick ref
    global_tropo = data.sel(lev=slice(700, 1100), lat=slice(-80, 80))
    # global_tropo = global_tropo.isel(time=-1)

    ####################################################################

    # Overview boxplots

    for subset in all_subsets():
        (lev, lev_slice), (lat, lat_slice) = subset
        print(lev, lat)

        data_subset = data.sel(lev=lev_slice, lat=lat_slice)
        # data_subset = data_subset.isel(time=-1)

        id_str = "{}_{}".format(lev, lat)

        for mode in args.mode:
            print(mode)
            fig = aerosol_boxplots(mode, data_subset)
            save_figure(_out_path("{}_{}".format(mode, id_str)),
                        fig=fig, qual='quick')
            plt.close(fig)

    n_pcts = mode_percentiles(global_tropo, 'n')
    mu_pcts = mode_percentiles(global_tropo, 'mu')

    print(n_pcts.swapaxes(0, 1)[['5', '95']])
    print(mu_pcts.swapaxes(0, 1)[['5', '95']])

    ####################################################################

    # Size distribution percentages - How many below a certain number
    # concentration?

    ns = [1e-4, 1e-3, 1e-2, 1e-1, 1.]
    data_set = global_tropo

    for n in ns:
        print("% of size distributions where N < {} for each mode".format(n))
        for mode in args.mode:
            print("{:10s}".format(mode), end=" - ")
            aer_data = data_set['n' + mode].values.ravel()

            pctile = percentileofscore(aer_data, n)
            print(pctile)
        print()

    ####################################################################

    # Comparison plots
    print("Comparison plots")
    for (mode, param) in product(args.mode, ['mu', 'n']):
        print("{} {}".format(mode, param))
        if ((mode.startswith("DST") or mode.startswith("SSLT")) and
           (param == 'mu')):
            continue
        g = cb(mode, param, all_subsets, data)
        save_figure(_out_path("{}_{}_comparison".format(mode, param)),
                    fig=g.fig, qual='quick')
        plt.close(g.fig)

    ocean_mask = extract_feature(data.isel(time=0, lev=0)['T'])
    data['ocean'] = (('lat', 'lon', ), ocean_mask)

    print("Comparison plots - maritime vs land")
    for (mode, param) in product(args.mode, ['mu', 'n']):
        print("{} {}".format(mode, param))
        if ((mode.startswith("DST") or mode.startswith("SSLT")) and
           (param == 'mu')):
            continue
        g = cmc(mode, param, data)
        save_figure(_out_path("{}_{}_ocean_vs_land_comparison"
                              .format(mode, param)),
                    fig=g.fig, qual='quick')
        plt.close(g.fig)

