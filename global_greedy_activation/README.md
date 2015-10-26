# Global, Greedy Activation Calculation

What modes are the most important for controlling aerosol activation? The script *calc_global.py* in this repo computes data to try to answer that question.

To do this, *calc_global* reads in some netCDF output from a CESM-MARC simulation which has already been pre-processed to contain detailed information on aerosol size distributions. It then subsets that output to identify situations where multiple aerosol modes were present, and uses those situations to drive [parcel model] calculations (using the updraft speed, temperature, aerosol mixture, etc.)

For each set of thermodynamic conditions and aerosol mixtures, a sequence of calculations are performed. First, an activation calculation is performed using *only* a single mode - that is, for the 12-mode mixture (4mode MARC subset + dust + sea salt), only one mode is present in the mixture. The calculation which produces the lowest Smax is said to be the "dominant" mode - it controls the activation dynamics. We retain that mode, and proceed to the second iteration calculation, which is for two-mode mixtures. The calculations are repeated, but now with the mode from the first iteration and a second mode from the remaining set. The combination which produces the lowest Smax is retained. 

The process is repeated until we run out of modes, or the variable `MAX_ITERS` is hit. At each step of the iteration, we track the mode which we add to the mixture, and its impact on Smax and the number concentration activated. The results are saved in two different formats, as well as all the intermediate calculations.

# *calc_global.py*

Rather than use a command-line interface, this script requires manual replacement of variables (all contained before the function definitions in the script):

- `N_THRESH`. `MU_THRESH`, and `ACT_THRESH` - tresholds for subsetting the aerosol populations in the MARC-CESM output file being scanned. Refer to the *iteration_setup_analysis* notebook for more details
- `act` - the name of the simulation to analyze
- `PARALLEL` - boolean toggle for running simulations in parallel
- `MAX_ITER` - maximum number of modes to include for breaking to "all-mode" calculation
- `ONLY_DISTS` - only select the paramteer sets and save them; don't run parcel model simulations 

# *iteration_setup_analysis* notebook

A handy notebook for visualizing how the choice of thresholds influenes the number of simulations available to study.

## Outputs

All output is stored in timestamped save directories from this one, e.g. `save_YYYYMMDD-hhmmss`. Each simulation set is saved as pickled, numbered output with corresponding output logs. Additionally, compiled CSV and netCDF files are saved in the directory, too. 

For analysis on machines other than legion, it's best to copy one of these directories. `arg_comp_ensemble` or a similar folder exists on my MBP for archiving. There are a *lot* of output files since each process generates individual outputs and logs. However, they compress down to very small archives using the following scheme:

- **concat_output.tar.gz**: `single_timestep_ti_*.nc`
- **concat_output_csv.tar.gz**: `single_timestep_ti_*.csv`
- **output_logs.tar.gz**: `out_****_******.log`
- **simulation_output.tar.gz**: `out_****_******`

## Log

- 9/25/2015 
    - re-extracted aerosol distributions from branch simulation 
    - updated `calc_global.py` script for new model simulations and toolkits
    - global simulations - [save_20150924-215548/]()
- 10/26/2015
    + Updated everything on MBP to most recent version
    + Backed up results on storage HDD


[parcel model]: http://www.github.com/darothen/parcel_model
