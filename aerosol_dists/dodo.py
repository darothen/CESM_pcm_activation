"""
Plotting pipeline management, using doit. To run, execute

$ doit [-n <number_of_jobs>]

"""

from plot_funcs import all_modes

#: File containing aerosol size distribution parameters
AEROSOL_DIST_FILE = "new_aerosol_dists.new_aerosols.nc"

## Size distribution parameters table
def task_dist_tables():
    """ Create a summary table of the different percentiles for number,
    mean size parameters for the different aerosol modes. """
    yield {
        'name': "dist_tables",
        'file_dep': ["dist_tables.py", AEROSOL_DIST_FILE, ],
        'actions': [["python", "dist_tables.py"], ],
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