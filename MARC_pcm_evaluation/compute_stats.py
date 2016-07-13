#!/usr/bin/env bash
""" Recompute summary statistics for a given chaos expansion.

"""

import os
import pandas as pd
import numpy as np
import sklearn.metrics as skm

from utils import summary_stats

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("sample_data", type=str, metavar="sampling.csv",
                    help="File containing sampling data")
parser.add_argument("parameterizations", type=str, nargs="+",
                    help="List of parameterizations to process from the"
                         " sampling results")
parser.add_argument("-e", "--exp-name", type=str, default='',
                    help="Name to label experiment; else infer from datafile")
parser.add_argument("-o", "--output", type=str,
                    help="Name of output file; else, will write "
                         "{exp_name}_stats_processed.p")


def result_key(output, name, val=0):
    """ Quickly alias a key into the sampling results dataset. """

    base_key = "%s_%s" % (output, name)
    if val > 0:
        base_key += "_%d" % val
    return base_key


def compute_stats_vs_parcel(df, output, param, output_parcel=None,
                            power10=False):
    if output_parcel is None:
        output_parcel = output

    # Figure out which columns in df to pull
    key_param = output+"_"+param
    key_parcel = output_parcel+"_parcel"

    # Grab and clean data
    data_df = pd.DataFrame({key_param: df[key_param],
                            key_parcel: df[key_parcel]})
    # Note that the data is all log data, so we're masking truly outlier
    # points here
    data_df[data_df > 20.] = np.nan
    data_df[data_df < -20.] = np.nan
    data_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    data_df.dropna(inplace=True)

    results = {}
    results['log10'] = summary_stats(data_df[key_param],
                                     data_df[key_parcel])
    data_df = 10.**(data_df)
    results['normal'] = summary_stats(data_df[key_param],
                                      data_df[key_parcel])

    return results

if __name__ == "__main__":

    args = parser.parse_args()

    # Load data
    print("Loading {}".format(args.sample_data))
    all_sampling_results = pd.read_csv(args.sample_data, index_col=0)
    if not args.exp_name:
        exp_name = os.path.basename(args.sample_data.split("_sampling")[0])
    else:
        exp_name = args.exp_name
    print("Successfully processed exp {}".format(exp_name))

    stats_df_mi = pd.MultiIndex.from_product(
        [["Smax", "Nderiv_Neq", "Nderiv_Nkn"],
         ["log10", "normal"],
         ["rmse", "nrmse", "mae", "r2", "mre", "mre_std"]],
        names=["result", "scaling", "stat"]
    )

    all_stats = []
    # for exp_name, exp in experiments.iteritems():

    print("Reading parameterizations...")
    for p in args.parameterizations:
        print("   " + p)

        stats = {}

        # Smax
        stats['Smax'] = compute_stats_vs_parcel(
            all_sampling_results, "Smax", p
        )

        # Nderiv vs Neq
        stats['Nderiv_Neq'] = compute_stats_vs_parcel(
            all_sampling_results, "Nderiv", p, output_parcel="Neq"
        )

        # Nderiv vs Nkn
        stats['Nderiv_Nkn'] = compute_stats_vs_parcel(
            all_sampling_results, "Nderiv", p, output_parcel="Nkn"
        )

        stats_vals = []
        for key in stats_df_mi:
            result, scaling, stat = key
            stats_vals.append(stats[result][scaling][stat])
        stats_df = pd.DataFrame(stats_vals, index=stats_df_mi,
                                columns=[p, ], )
        all_stats.append(stats_df.T)
    # END INDENT

    all_df = pd.concat(all_stats)
    all_df.index.name = "scheme"

    # Write output
    if not args.output:
        out_fn = "{}_stats_processed.p".format(exp_name)
    else:
        out_fn = args.output
    all_df.to_pickle(out_fn)
