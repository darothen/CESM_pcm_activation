#!/usr/env/bin python
"""
Summarize MARC aerosol size distributions in text-tabular format.

"""

import pandas as pd
import xray

from plot_funcs import mode_percentiles

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("aerosol_ds", type=str, metavar="aerosol_data.nc",
                    help="Input file containing the aerosol size"
                         " distribution parameters")

if __name__ == "__main__":

    args = parser.parse_args

    # Read in aerosol dataset
    data = xray.open_dataset(args.aerosol_ds, decode_times=False)

    # Nudge times to the year 2000
    times = data.coords.to_dataset().time
    times += 2000.*365
    data = xray.decode_cf(data)

    # Global troposphere slice for quick ref
    global_tropo = data.sel(lev=slice(700, 1100), lat=slice(-80, 80))
    global_tropo = global_tropo.isel(time=slice(10, None))
    # global_tropo = global_tropo.isel(time=-1)

    # Compute percentiles
    n_pcts = mode_percentiles(global_tropo, 'n')
    mu_pcts = mode_percentiles(global_tropo, 'mu')

    # Output tables
    pd.set_option('display.max_columns', 9999)
    pd.set_option('precision', 3)

    with open("figs/aerosol_dists_number.table", 'w') as f:
        print(n_pcts.to_string(), file=f)
    with open("figs/aerosol_dists_mu.table", 'w') as f:
        print(mu_pcts.to_string(), file=f)