"""
Define some basic styling, context, and data contexts for the
plotting scripts
"""
import os, pickle, sys

from collections import namedtuple
from itertools import product, cycle

from string import ascii_lowercase

import seaborn as sns
sns.set(context='talk', style='ticks')

import pandas as pd
import numpy.random as npr

PAPER_DIR = "/Users/daniel/Documents/Writing/pce_single_mode"

## Choose a data source
DATA_LOG_DIR = os.path.join(PAPER_DIR, "data_20k_log")
DATA_LIN_DIR = os.path.join(PAPER_DIR, "data_10k_lin")

def get_stats(kind='lin', result='Smax', scaling='normal', override=False):
    """ Retrieve the statistics from the sampling experiments. 

    Parameters
    ----------
    kind : str
        'lin', 'log'
    result : str
        'Smax', 'Neq', 'Nderiv_Neq', etc.
    scaling : str
        'normal' or 'log10'

    """

    stats_fn = "SM_stats.p"
    if kind == 'log':
        with open(os.path.join(DATA_LOG_DIR, stats_fn), 'r') as f:
            stats = pickle.load(f)
    else:
        with open(os.path.join(DATA_LIN_DIR, stats_fn), 'r') as f:
            stats = pickle.load(f)

    if override: return stats

    s = stats.copy()
    # Choose variable
    s = s.T.xs(result, level='result').T
    # Choose sampling
    s = s.T.xs(scaling, level='scaling').T

    return s

def get_data(kind='lin', which='tidy'):
    """ Retrieve the data from the sampling experiments. 

    Hard-coded to assume that the 'log' experiment has twice as many
    samples as the 'linear' one.

    Parameters
    ----------
    kind : str
        'lin', 'log', or 'mix'
    which : str
        'tidy', 'results'

    """
    data_fn = "SM_sampling_%s.csv" % which

    if kind in ['log', 'both']:
        data_log = pd.read_csv(os.path.join(DATA_LOG_DIR, data_fn),
                               index_col=0)
        data_log = data_log.loc[npr.choice(data_log.index, len(data_log)/2, 
                                           replace=False)]
        data_log = data_log.reset_index(drop=True)
    if kind in ['lin', 'both']:
        data_lin = pd.read_csv(os.path.join(DATA_LIN_DIR, data_fn),
                               index_col=0)
        data = data_lin

    if kind == 'both':
        data = pd.concat([data_log, data_lin], ignore_index=True)
    else: 
        data = data_log if kind == 'log' else data_lin

    return data

EXP_PREFIX = 'SM'
EXP_NAMES = ['OLS', 'LARS', 'LASSO']
EXP_ORDERS = [2, 3, 4, 5]
PARAM_NAMES = ['ARG', 'MBN']

def get_compare_colors():
    return cycle(sns.color_palette('colorblind', 5))
hue_order = EXP_NAMES + PARAM_NAMES

## Alias for labeling subplots
def get_letters():
    return cycle(ascii_lowercase)
letter_coords = 0.05, 0.9
letter_fontdict = dict(fontsize=14, fontweight='bold')



## Long name mapping
def fix_varname(varname):
    if "_" in varname:
        pcm, lev = varname.strip().split("_")
        return pcm
    elif varname in ['ARG', 'MBN', 'parcel']:
        return {"ARG": "ARG2011", "parcel": "Numerical", "MBN": "MBN2014"}[varname]
    elif varname == ['MARC_ols']:
        return "OLS"
    else:
        return varname

short_parameter_names = {
    "logV": "$\log_{10} V$",
    "P": "$P$",
    "T": "$T$",
    "logmu": "$\log_{10} \mu$",
    "logN": "$\log_{10} N$",
    "sigma": "$\sigma$",
    "kappa": "$\kappa$",
    "accom": "$a_c$"
}

long_parameter_names = {
    "logV": "log10(Updraft Speed [m/s])",
    "P": "Pressure (Pa)",
    "T": "Temperature (K)",
    "logmu": "log10(Median Radius [micron])",
    "logN": "log10(Number Concentration [1/cm$^3$])",
    "sigma": "Geometric Std Dev",
    "kappa": "Hygroscopicity",
    "accom": "Condensation Coeficient"
}

## Functions
fn_fix = lambda s: s.replace("_", "\_")

