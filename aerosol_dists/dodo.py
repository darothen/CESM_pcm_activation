"""
Plotting pipeline management, using doit. To run, execute

$ doit [-n <number_of_jobs>]

"""

import sys
import warnings

from plot_funcs import all_modes
try:
    import marc_analysis
except ImportError:
    warnings.warn("marc_analysis package not found! ")
    sys.exit()

#: Original instantaneous aerosol output
MARC_AEROSOL_OUT = "/Volumes/legion_home/CESM_archive/act_aie_aci/aerosols/aerosol_dists_PD.sample.nc"

#: File containing aerosol size distribution parameters
AEROSOL_DIST_FILE = "aerosol_dists.nc"

## Compute aerosol size distribution data
def task_calc_aerosol_dists():
    """ Compute the aerosol size distribution parameters from the CESM/MARC
    output tapes.

    Note
    ----
    This is a big calculation if we try to stream the data over an internet
    connection. It's easiest to log into legion and simply execute

        darothon@legion$ calc_aerosol $MARC_AEROSOL_OUT $AEROSOL_DIST_FILE

    and then copying the results here. Then, this task can be ignored through
    doit

        daniel@roth-mbp$ doit ignore calc_aerosol_dists

    """
    warnings.warn("STOP! run this command manually and then mark as complete")
    # yield {
    #     'name': "calc_aerosol_dists",
    #     'targets': [AEROSOL_DIST_FILE, ],
    #     'file_dep': [MARC_AEROSOL_OUT, ],
    #     'actions': [["calc_aerosol", MARC_AEROSOL_OUT, AEROSOL_DIST_FILE, ]],
    # }

## Size distribution parameters table
def task_dist_tables():
    """ Create a summary table of the different percentiles for number,
    mean size parameters for the different aerosol modes. """
    yield {
        'name': "dist_tables",
        'file_dep': ["dist_tables.py", AEROSOL_DIST_FILE, ],
        'actions': [["python", "dist_tables.py", AEROSOL_DIST_FILE, ], ],
    }

## Comparison plots
def task_compare_mode():
    """ Generate boxplots of the aerosol size distribution parameters
    for different latitude/level slices, etc. """

    for mode in all_modes:
        task_dict = {
            'name': mode,
            'file_dep': ['compare_dists.py', AEROSOL_DIST_FILE, ],
            'actions': [['python', 'compare_dists.py', AEROSOL_DIST_FILE,
                         '-m', mode], ],
        }
        yield task_dict