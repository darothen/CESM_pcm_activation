#!/usr/bin/env python
""" Greedy aerosol acitvation calculations

Restrict evaluation to grid cells which meet the following criteria,
on top of the

1. Grid levels 23-26, or about 850 - 950 mb
    - taken care of in extraction notebook
2. Mid-latitudes, -60 < lat < 60
3. Activation *actually occurred* in that timestep
    - sum all the {MODE}ACT values, and ensure nonzero

As a test, try using xray for the extraction; we'll have to

TODO: add restart functionality, where user can specify SAVE_DIR and
      the script will lookup succesful simulations from cache

"""

from __future__ import print_function

import os
import subprocess
import time

import parcel_model

import numpy
import numpy as np
import pandas as pd
import xray

mode_dict = {
    'OC':  { 'sigma': 2.0, 'rho': 2000., 'kappa': 1e-10 },
    'MOS': { 'sigma': 2.0, },
    'BC':  { 'sigma': 2.0, 'rho': 2000., 'kappa': 1e-10 },
    'MBS': { 'sigma': 2.0, 'rho': 2000., 'kappa': 0.507 },
    'NUC': { 'sigma': 1.59, 'rho': 1800., 'kappa': 0.507 },
    'AIT': { 'sigma': 1.59, 'rho': 1800., 'kappa': 0.507 },
    'ACC': { 'sigma': 1.59, 'rho': 1800., 'kappa': 0.507 },
    'DST': { 'mus': [0.16, 0.406, 0.867, 1.656], 'kappa': 0.14,
             'sigmas': [1.4,]*4 },
    'SSLT':{ 'mus': [0.50, 2.000, 5.000,15.000], 'kappa': 1.16 ,
             'sigmas': [1.59, 1.37, 1.41, 1.22] },
}
modes = mode_dict.keys()
rho_OC = mode_dict['OC']['rho']
rho_SO4 = mode_dict['ACC']['rho']
kappa_OC = mode_dict['OC']['kappa']
kappa_SO4 = mode_dict['ACC']['kappa']

MARC_modes = ['AIT', 'ACC', 'MOS', 'MBS']
DST_modes  = ["DST%02d" % d for d in range(1, 5)]
SSLT_modes = ["SSLT%02d" % d for d in range(1, 5)]
all_modes  = MARC_modes + DST_modes + SSLT_modes

#: Experiment name
NAME = "single_timestep"

#: Diagnosed in iteration_setup_analysis.ipynb
N_THRESH, MU_THRESH = 1.0, 2e-3

#: Aerosol instantaneous samples
INPUT_FN = "new_aerosol_dists.new_aerosols.nc"

#: Random/parcel model
RANDOM = False

#: Calculations in parallel?
PARALLEL = True

#: Number of iterations before final, all aerosol calculation
MAX_ITER = 4

#: Only output satisfactory distributions from the input dataset
ONLY_DISTS = False

#: Dataset subsampling
LEVEL_SUB = slice(500, 1100)
LAT_SUB = slice(-70, 70)
TIME_ISUB = slice(-10, None)

def get_git_versioning():
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

def smax_nact_calc(aerosol_modes, V, T, P):
    import parcel_model as pm

    if RANDOM:
        print("<<< RANDOM >>>")
        Smax = numpy.random.random()
        eq = numpy.random.random()
        kn = numpy.random.random()
        return Smax, eq, kn

    output_dt = 0.1
    solver_dt = 10.0
    t_end = 1800.
    accom = 0.1

    try:
        model = pm.ParcelModel(aerosol_modes, V, T, -0.0, P, accom=accom,
                               truncate_aerosols=True, console=False)
        p, a = model.run(t_end=t_end, output_dt=output_dt, solver_dt=solver_dt,
                         max_steps=2000, solver='cvode', output='dataframes',
                         terminate=True, terminate_depth=10.0)
        Smax = p['S'].max()
        T_fin = p['T'].iloc[-1]
        parcel_sucess = True

        eq, kn = 0., 0.
        for aer in aerosol_modes:
            rs = a[aer.species].iloc[-1].values
            af_eq, af_kn, _, _ = pm.binned_activation(Smax, T_fin,
                                                      rs, aer)
            # Use binned totals for consistency
            eq += af_eq*numpy.sum(aer.Nis*1e-6)
            kn += af_kn*numpy.sum(aer.Nis*1e-6)

    except parcel_model.ParcelModelError:
        Smax, eq, kn = 0., 0., 0.
        parcel_sucess = False

    return Smax, eq, kn

def exec_run(args):
    from timeit import default_timer as timer
    import parcel_model as pm

    ti, i, ind, row = args

    # output file stream
    f = open(os.path.join(SAVE_DIR, "out_%04d_%06d.log" % (ti, ind)), 'w')
    if PARALLEL: # print to file
        def ptf(*s):
            print(*s, file=f); f.flush()
    else:
        def ptf(*s):
            print(*s); f.flush()

    V, T, P,                            \
    muAIT, muACC, muMOS, muMBS,         \
    nAIT, nACC, nMOS, nMBS,             \
    nDST01, nDST02, nDST03, nDST04,     \
    nSSLT01, nSSLT02, nSSLT03, nSSLT04, \
    kappaMOS, time, lev, lat, lon       = row

    ptf( "\n   Simulation (%d - %d/%d)" % (ti, i, ind),  )
    ptf( "   -------------------------------", )
    ptf( "       ACC:", muACC, nACC )
    ptf( "       AIT:", muAIT, nAIT )
    ptf( "       MBS:", muMBS, nMBS)
    ptf( "       MOS:", muMOS, nMOS, kappaMOS)
    ptf( "      nDST:", nDST01, nDST02, nDST03, nDST04)
    ptf( "     nSSLT:", nSSLT01, nSSLT02, nSSLT03, nSSLT04)
    ptf( "    (time, lev, lat, lon) = (%d, %d, %d, %d)" % (time, lev, lat, lon))

    # Setup the aerosol
    aerosol_modes = []
    aer_dict = {}

    aer_dict["ACC"] = pm.AerosolSpecies("ACC",
                            pm.Lognorm(mu=muACC, sigma=mode_dict['ACC']['sigma'],
                                       N=nACC),
                            kappa=mode_dict['ACC']['kappa'], bins=100)
    aerosol_modes.append(aer_dict["ACC"])

    aer_dict["AIT"] = pm.AerosolSpecies("AIT",
                            pm.Lognorm(mu=muAIT, sigma=mode_dict['AIT']['sigma'],
                                       N=nAIT),
                            kappa=mode_dict['AIT']['kappa'], bins=100)
    aerosol_modes.append(aer_dict["AIT"])

    aer_dict["MBS"] = pm.AerosolSpecies("MBS",
                            pm.Lognorm(mu=muMBS, sigma=mode_dict['MBS']['sigma'],
                                       N=nMBS),
                            kappa=mode_dict['MBS']['kappa'], bins=100)
    aerosol_modes.append(aer_dict["MBS"])

    aer_dict["MOS"] = pm.AerosolSpecies("MOS",
                            pm.Lognorm(mu=muMOS, sigma=mode_dict['MOS']['sigma'],
                                       N=nMOS),
                            kappa=kappaMOS, bins=100)
    aerosol_modes.append(aer_dict["MOS"])

    for i, (dst, sslt) in enumerate(zip([nDST01, nDST02, nDST03, nDST04],
                                        [nSSLT01, nSSLT02, nSSLT03, nSSLT04])):
        n = "%02d" % (i+1, )
        for name, dat in zip(['DST', 'SSLT'], [dst, sslt]):
            aer_dict[name+n] = pm.AerosolSpecies(name+n,
                                    pm.Lognorm(mu=mode_dict[name]['mus'][0],
                                               sigma=mode_dict[name]['sigmas'][0],
                                               N=dat),
                                    kappa=mode_dict[name]['kappa'], bins=50)
            aerosol_modes.append(aer_dict[name+n])

    it = 0
    modes  = [] # name of mode added at iteration i
    smaxes = [] # smax achieved by mixture at iteration i
    nacts  = [] # total activated at iteration i
    active_modes = [] # register of which modes are present in aerosol mixture

    # timer setup
    start = timer()

    # greedy iteration - begin
    ptf("\n   Begin greedy iterations")
    ptf("   -----------------------")

    while it < MAX_ITER:

        iter_start = timer()

        ptf("      %02d) Active modes" % (it+1, ))

        min_smax = 1e20
        min_mode = None
        min_index = None
        for mi, mode in enumerate(aerosol_modes):
            current_modes = active_modes + [mode, ]

            smax, eq, kn = smax_nact_calc(current_modes, V, T, P)

            # Status update on coneole
            ptf("          (%s) - %3.2e" %
                (",".join([m.species for m in current_modes]), smax))

            # Handle error conditions:
            # 1) Model failed; Smax <= 0. --> Continue to next
            if smax <= 0.:
                continue

            # else, we've succeeded!
            if smax < min_smax:
                min_smax = smax
                min_mode = mode
                min_nact = [eq, kn]
                min_index = mi

        dt = timer() - iter_start
        ptf("      Elapsed time: %5.1f s" % dt)
        ptf("      ----> %s (smax = %3.2e)" % (min_mode.species, min_smax))
        ptf("      -----------------------")

        modes.append(min_mode.species)
        smaxes.append(min_smax)
        nacts.append(min_nact)
        active_modes.append(min_mode)

        del aerosol_modes[min_index]

        it += 1

    # Finish with all the modes
    ptf( "      Final simulation with all modes: ")
    current_modes = active_modes + aerosol_modes
    smax, eq, kn = smax_nact_calc(current_modes, V, T, P)
    ptf( "             smax = %3.2e" % smax )
    modes.append('all')
    smaxes.append(smax)
    nacts.append([eq, kn])

    # Convert nacts to ndarray; it's not useful right now
    nacts = numpy.array(nacts)
    eqs = nacts[:, 0]
    kns = nacts[:, 1]

    # Save an intermediate file, as a cache
    with open(os.path.join(SAVE_DIR, "out_%04d_%06d" % (ti, ind)), 'wb') as f_save:
        numpy.save(f_save, [modes, smaxes, eqs, kns])

    dt_tot = timer() - start
    ptf("   Total time: %3.1f minutes" % (dt/60., ))

    # close the log
    f.close()

    # Return
    return [modes, smaxes, eqs, kns]

if __name__ == "__main__":

    print("CALC_GLOBAL: greedy single-mode activation calculations")
    input_fn = INPUT_FN
    print("    analyzing file %s" % input_fn)

    # Construct the save directory
    timestamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    SAVE_DIR = os.path.join(os.getcwd(), "save_%s" % timestamp)
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR)
    print("    saving into directory %s" % SAVE_DIR)

    if PARALLEL:
        from ipyparallel import Client

        # PARALLEL PROCESSING CONTROL
        client = Client(profile="legion")
        dv = client[:]

        with dv.sync_imports():
            import os
            import parcel_model
            import numpy

        dv['exec_run'] = exec_run
        dv['smax_nact_calc'] = smax_nact_calc

        dv['mode_dict'] = mode_dict
        # Correct save_dir for network mount
        SAVE_DIR = os.path.join("/net/legion", SAVE_DIR)
        dv['SAVE_DIR'] = SAVE_DIR
        dv['MAX_ITER'] = MAX_ITER
        dv['PARALLEL'] = PARALLEL
        dv['RANDOM'] = RANDOM

        view = client.load_balanced_view()
        print("   parcel model calculations will be done in parallel")

    data = xray.open_dataset(input_fn,
                             decode_times=False # we don't care about datetimes
    )
    # Subset data fields
    data = data.sel(lev=LEVEL_SUB, lat=LAT_SUB)
    data = data.isel(time=TIME_ISUB)

    # Alias for extracting data underlying an xray data set
    # given its key
    def d(field, ravel=False):
        data_array = data[field].data
        if ravel:
            return data_array.ravel()
        return data_array

    # Build a multi-dimensional index set for extracting
    time_4d, lev_4d, lat_4d, lon_4d = \
        np.meshgrid(*map(d, ['time', 'lev', 'lat', 'lon']),
                    indexing='ij')

    print("    processing some useful values...")

    lt = lambda a, b: a & b
    tot_samples = 0

    # Iterate over time slices to chunk the calculations
    print("\nIterating over timeslices...")
    print("----------------------------")
    for ti in range(len(d('time'))):
        t = d('time')[ti]

        print("   {:03d}): {:8.1f}".format(ti, t))

        # mapping =  [ ( aeract[ti], ACT_THRESH), ]
        mapping = []
        mapping += [ ( d('n%s' % mode)[ti], N_THRESH ) for mode
                                                       in MARC_modes ]
        mapping += [ ( d('mu%s' % mode)[ti], MU_THRESH ) for mode
                                                        in MARC_modes ]

        fields, threshholds = zip(*mapping)
        mask = reduce(lt, [(field > thresh) for field, thresh
                                            in zip(fields, threshholds)])
        print("   Number of valid samples: ", np.count_nonzero(mask))

        # Now, the data to use can be indexed by simply calling
        # d(field)[ti, mask]

        # Gather the values to use as parameter sets. This isn't elegant, but
        # it sure as hell gets the job done straightaway.
        # `field_names` is the order of the arguments to pass to each iterative
        # solver
        field_names = ["WSUB", "T", "P", ] + \
                      ["mu%s" % mode for mode in MARC_modes] + \
                      ["n%s" % mode for mode in all_modes] + \
                      ["kappaMOS", ]
        rows = np.array([d(field)[ti, mask] for field in field_names] +
                        [time_4d[ti, mask], lev_4d[ti, mask], lat_4d[ti, mask],
                         lon_4d[ti, mask]] )
        df = pd.DataFrame(rows.T,
                          columns=field_names + ["time", "lev", "lat", "lon"]).dropna()
        assert len(df) == np.count_nonzero(mask)

        # Short circuit to just get dists
        if ONLY_DISTS:
            df.to_csv(os.path.join(SAVE_DIR, "./%s_ti_%04d.csv" % (NAME, ti) ))
            continue

        parallel_rows = [[ti, i, ind, row.values] for i, [ind, row] in enumerate(df.iterrows())]
        if PARALLEL:
            results = view.map_async(exec_run, parallel_rows, ordered=True)
            results.wait_interactive()

        else:
            results = []
            for row in parallel_rows:
                results.append(exec_run(row))

        # TODO: I break down the results into their components and manually
        # construct the reuslts dataframe which I save to disk. There must be
        # a better way to construct the results array, using mixed data types,
        # but each result is its own matrix with mixed datatypes!

        mode_results = np.array([r[0] for r in results])
        numerical_results = np.array([r[1:] for r in results])

        # Construct a dataframe to store the results in "wide" form
        mode_cols = [ "MODE_%02d" % (i+1, ) for i in range(MAX_ITER) ] + ["MODE_all", ]
        smax_cols = [ "SMAX_%02d" % (i+1, ) for i in range(MAX_ITER) ] + ["SMAX_all", ]
        eq_cols   = [ "EQ_%02d" % (i+1, ) for i in range(MAX_ITER) ]   + ["EQ_all", ]
        kn_cols   = [ "KN_%02d" % (i+1, ) for i in range(MAX_ITER) ]   + ["KN_all", ]

        results_df = pd.DataFrame(columns=mode_cols+smax_cols+eq_cols+kn_cols,
                                  index=df.index)
        results_df[mode_cols] = mode_results
        results_df[smax_cols] = numerical_results[:,0,:]
        results_df[eq_cols]   = numerical_results[:,1,:]
        results_df[kn_cols]   = numerical_results[:,2,:]

        combined_df = df.join(results_df).dropna()
        combined_df.to_csv(os.path.join(SAVE_DIR,
                                        "./%s_ti_%04d.csv" % (NAME, ti) ))

        # Save in a tidy netCDF file via xray
        nrecs, niters = mode_results.shape
        ds = xray.Dataset(
            { 'MODE': (('index', 'iter'), mode_results),
              'SMAX': (('index', 'iter'), numerical_results[:,0,:]),
              'EQ': (('index', 'iter'), numerical_results[:,1,:]),
              'KN': (('index', 'iter'), numerical_results[:,2,:]),   },
            coords = {
                'iter': range(1, niters+1), 'index': range(nrecs),
            },
        )
        ds.attrs['input_file'] = input_fn
        ds.attrs['timestamp'] = timestamp
        ds.attrs['ti'] = ti
        ds.attrs['n_thresh'] = N_THRESH
        ds.attrs['mu_thresh'] = MU_THRESH
        # ds.attrs['git_commit'] = get_git_versioning()
        ds.attrs['numpy_version'] = np.__version__
        ds.attrs['pandas_version'] = pd.version.version
        ds.attrs['xray_version'] = xray.version.version

        params_ds = xray.Dataset.from_dataframe(df)
        combined_ds = xray.auto_combine([ds, params_ds])

        # Increment index so order is sequential over all chunks
        ds['index'] = ds.index + tot_samples

        # Save to disk
        combined_ds.to_netcdf(os.path.join(SAVE_DIR, "%s_ti_%04d.nc" % (NAME, ti)))

        # Record of how many samples we've processed
        tot_samples += len(df.index)

    data.close()


