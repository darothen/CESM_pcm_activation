"""
Utilities for visualizing output from the greedy activation
simulations.
"""
from __future__ import print_function

import os, glob

import numpy as np
import pandas as pd
import xray

import geopandas
from shapely.geometry import asPoint, MultiPolygon
from shapely.prepared import prep
oceans = geopandas.read_file("ne_110m_ocean.shp")
oceans = MultiPolygon(oceans.geometry.values.tolist())
oceans = prep(oceans)

def shift_lons(lons):
    mask = lons > 180
    lons[mask] = -(360. - lons[mask])
    return lons

def _is_in_ocean(p, oceans=oceans):
    return oceans.contains(p)

def calc_ocean_mask(dataset, oceans=oceans, pt_return=False,
                    longitude='lon', latitude='lat'):

    lons, lats = dataset[longitude], dataset[latitude]
    if isinstance(dataset, (xray.Dataset, xray.DataArray)):
        lons.values = shift_lons(lons.values)
    else:
        lons = shift_lons(lons)

    points = [asPoint(point) for point in np.column_stack([lons, lats])]
    if pt_return:
        return points

    in_ocean = [_is_in_ocean(p, oceans) for p in points]
    return in_ocean

def load_experiment(exp_dir, name="single_timestep", format='csv'):
    """ Load the results from a complete experiment into a DataFrame.

    Parameters
    ----------
    exp_dir : str
        The path to the directory containing the output CSV files
        from the simulations.
    name : str
        The name of the experiment files; default is "single_timestep"
    format : str
        Either "csv" or "nc" for loading the correct input

    """

    fns = sorted(glob.glob(os.path.join(exp_dir, "%s*.%s" % (name, format))))
    print("Found %d files" % len(fns))

    if format == 'csv':

        dfs = []
        print("Reading...")
        for fn in fns:
            print ("   ", fn)
            dfs.append(pd.read_csv(fn, index_col=0))

        df = pd.concat(dfs, ignore_index=True)

        return df

    elif format == 'nc':

        print("Reading...")
        ds = xray.open_mfdataset(fns)
        ds.set_coords(["lat", "lon", "lev"], inplace=True)
        return ds

    else:
        raise ValueError("Format should either be 'nc' or 'csv'.")