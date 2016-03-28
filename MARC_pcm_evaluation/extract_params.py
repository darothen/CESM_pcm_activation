#!/usr/bin/env python
"""
Extract parameters from a CESM/MARC simulation to use to compare a chaos
expansion to parcel model output.

Notes
-----


"""
from __future__ import print_function

import os, pickle, warnings
from itertools import product

import numpy as np
import pandas as pd

import xray

import argparse
parser = argparse.ArgumentParser(description="Extract aerosol activation calculation parameters" +
                                             " from a CESM/MARC simulation")
parser.add_argument("exp_name", type=str, help="name of PCM experiment to reference")
parser.add_argument("-s", "--subset", action='store_true',
                    help="Subset CESM data based on PCM bounds")
parser.add_argument("-f", "--cesm_file", type=str, required=True,
                    help="CESM/MARC output data to extract from")
parser.add_argument("-n", "--n", type=int, default=10000,
                    help="Number of samples to extract from the dataset")
parser.add_argument("-o", "--output", type=str, default="CESM_extract.csv",
                    help="Name of output file")

#: Aerosol modes to extract
MARC_modes = ["ACC", "MOS", "MBS", ]
DST_modes = ["DST01", "DST02", ]
SSLT_modes = ["SSLT01", ]

#: Meteorology fields to extract
meteo_fields = ["Q", "P", "T", "WSUB", "ORO"]

#: water vapor threshold
QCSMALL = 1e-5

#: Mapping of PCM var shorthand to CESM fields
mapping = {
    "WSUB": "V",
}

if __name__ == "__main__":

    args = parser.parse_args()
    print("Referencing PCM %s" % args.exp_name)
    print("Extracting (%d) samples from %s" % (args.n, args.cesm_file))
    print("--"*20 + "\n")

    ###########################################################

    ## Load in the experiment dictionary, which contains information on
    ## the parameter space structure and contents
    exp_dict = pickle.load(open("%s_exp.dict" % args.exp_name, 'r'))

    directive_base = exp_dict['directive_base']
    variables = exp_dict['variables']
    responses = exp_dict['responses']

    print("Reading in activation parameter fields")

    ds = xray.open_dataset(args.cesm_file, decode_cf=False)

    # Subset the original dataframe
    ds_subset = ds[
        ["n"+aer for aer in MARC_modes+DST_modes+SSLT_modes] +
        ["mu"+aer for aer in MARC_modes] +
        ["kappaMOS", ] +
        meteo_fields
    ]

    # Further subset by narrowing the lev and lat bands
    ds_subset = ds_subset.sel(lev=slice(700., 1000.),
                              lat=slice(-60., 60.))

    # Re-name necessary fields
    ds_subset.rename(mapping, inplace=True)

    ## Subset the variables to match the PCM setup
    print("Mapping the PCM vars to the CESM data")
    pcm_vars = [v[0] for v in variables]
    for var in pcm_vars:
        logvar = "log" in var
        print("   ", var, " -> ", end=' ')

        ## Some logic to map the PCM varnaming scheme to the CESM modes
        if len(var) == 1:
            var_short = var
        else:
            bits = var.split("_")
            if len(bits) == 1 and logvar:
                var_short = bits[0][3:]
            elif len(bits) > 1 and logvar:
                prefix, var_short = var.split("_")
                prefix = prefix[3:]
                var_short = prefix.lower()+var_short
            else:
                var_short = var.replace("_", "")
        if var_short in mapping:
            var_short = mapping[var_short]

        if logvar:
            ds_subset[var_short] = np.log10(ds_subset[var_short])
        print(var_short)
        ds_subset.rename({var_short: var}, inplace=True)

    # Mask where QC is too small, kappa is invalid
    ds_masked = ds_subset.where(ds_subset['Q'] > QCSMALL)

    # Convert to dataframe, drop NA
    print("Converting to DataFrame...")
    df = ds_masked.to_dataframe().reset_index().dropna()
    n_tot = len(df)

    row = "   {:<14s} {:<12d}"
    print(row.format("valid data points", n_tot))
    if args.subset:
        warnings.warn("NOT IMPLEMENTED - analyze w/o bounds beforehand")
        # print("Subsetting data based on PCM bounds")
        # for v in variables:
        #     pcm_var = v[0]
        #     _, lo, hi = v[3]
        #     mask = mask & (data[data_mapping[pcm_var]] > lo) & \
        #                   (data[data_mapping[pcm_var]] < hi)
        #     print(row.format(pcm_var, len(data['Q'][mask])))
        # n_tot = len(data['Q'][mask])

    ## Sample the data
    print("Sampling the data")
    random_df = df.sample(args.n)

    ## Finally print a table describing the ranges of the extracted vars
    print("\n" + "--"*30)
    print("Distributions of extracted data" + "\n")
    header = " "*14 + " {:<10s}"*4
    print(header.format('min', 'median', 'mean', 'max'))

    row = "{:<14s}" + " {:<10.2f}"*4
    for var in random_df:
        var_dat = random_df[var]
        try:
            row_prnt = [var, np.min(var_dat), np.median(var_dat),
                        np.mean(var_dat), np.max(var_dat)]
            print(row.format(*row_prnt))
        except:
            print("skipping " + var)
    print()

    print("Saving to file %s" % args.output)
    random_df.to_csv(args.output)