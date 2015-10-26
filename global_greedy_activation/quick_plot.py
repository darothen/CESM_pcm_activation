#!/usr/bin/env python
""" 
"""

from __future__ import print_function

import argparse
import os, sys, time
sys.path.insert(0, "../")

from case_setup import *

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

#: Diagnosed in iteration_setup_analysis.ipynb
# N_THRESH, MU_THRESH = 1., 2e-3 
N_THRESH, MU_THRESH = 1., 2e-3 
ACT_THRESH = 100.

#: Activation case
act = "pcm_min_smax"

class QuickPlot(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description="Quickly plot iteration sample info",
            usage="quick_plot <experiment> <command> [<args>]"
        )

        parser.add_argument("experiment", help="Experiment to use for data")
        parser.add_argument("command", help="Subcommand to run")
        parser.add_argument("-n","--N", type=float, default=N_THRESH, 
                            help="N threshold")
        parser.add_argument("-m","--mu", type=float, default=MU_THRESH, 
                            help="mu threshold")
        parser.add_argument("-a","--act", type=float, default=ACT_THRESH,
                            help="act threshold")
        args = parser.parse_args(sys.argv[1:])

        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            sys.exit(1)

        self._args = args
        self._init_data(args.experiment) # set -> self._df

        # dispatch pattern -> invoke method with arg name
        getattr(self, args.command)()

    def test(self):
        print("This is just a test")

    def _init_data(self, act):
        input_fn = "data/%s.aerosols.nc" % act
        print("    analyzing file %s" % input_fn )

         # Construct the save directory
        timestamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())

        data = xray.open_dataset("data/%s.aerosols.nc" % act,
                                 decode_times=False # we don't care about datetimes
        )
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
        dt = 1800. # seconds
        aeract = -1.*d('nAERACT')*d('RHO')*dt*1e-6

        lt = lambda a, b: a & b
        tot_samples = 0

        # Iterate over time slices to chunk the calculations
        print("\nIterating over timeslices...")
        print("----------------------------")
        dfs = []
        for ti in xrange(len(d('time'))):
            t = d('time')[ti]

            print("   {:03d}): {:8.1f}".format(ti, t))

            mapping =  [ ( aeract[ti], self._args.act), ]
            mapping += [ ( d('n%s' % mode)[ti], self._args.N ) for mode
                                                           in MARC_modes ]
            mapping += [ ( d('mu%s' % mode)[ti], self._args.mu ) for mode 
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

            dfs.append(df)

        data.close()

        self._df = pd.concat(dfs, ignore_index=True)
        pf_args = dict(n_samples=len(self._df), 
                       N=self._args.N, mu=self._args.mu, act=self._args.act)
        print(("\n...found {n_samples:d} samples for thresholds" + 
               " N={N:3.1e}, mu={mu:3.1e}, act={act:3.1e}\n").format(**pf_args))

    def dists(self):
        import matplotlib.pyplot as plt
        import seaborn as sns

        print("Plotting distributions for all parameters...")

        keys = [k for k in self._df.keys() if not "_" in k]
        n_plots = len(keys)
        n_cols = 4
        n_rows = n_plots / n_cols

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*3.5, n_rows*2))
        for i, (key, ax) in enumerate(zip(keys, axs.ravel())):
            
            kde_args = dict(ax=ax, shade=True)
            if i > 0: kde_args['legend'] = False
            
            data = self._df[key]
            if key.startswith("n") or key.startswith("mu"):
                data = np.log10(data)
                key = "log10(%s)" % key
            
            sns.kdeplot(data, color='k', **kde_args)
            ax.set_xlabel(key)
            plt.setp(ax.get_xticklabels(), rotation=20)
            
            if (i % n_cols) == 0:
                ax.set_ylabel("density function")

        sns.despine(fig=fig)
        plt.tight_layout()

        plt.show()

    def map(self):
        import matplotlib.pyplot as plt
        import seaborn as sns

        import cartopy.crs as ccrs
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


        print("Plotting global distribution of samples...")

        dlon, dlat = 2.5, 1.9
        npoints = len(self._df)
        df_map = self._df.copy()

        def rand_fact(npoints, width):
            """ Rescale a number from [0, 1) to [-width, width) """
            points = 2.*np.random.random(npoints) - 1. 
            points *= width
            return points

        l1 = df_map['lon'].copy()
        l2 = df_map['lat'].copy()

        df_map['lon'] = df_map['lon'] + rand_fact(npoints, dlon/2.)
        df_map['lat'] = df_map['lat'] + rand_fact(npoints, dlat/2.)

        lon, lat = df_map.lon, df_map.lat
        # Correct lon: 
        # 1) some values may be < 0 or > 360, map these into [0, 360]
        lon[lon < 0] = lon[lon < 0] + 360.
        lon[lon > 360] = lon[lon > 360] - 360.
        # 2) map from [0, 360] -> [-180, 180]
        lon -= 180.

        df_map['lon'] = lon[:]     

        proj = ccrs.PlateCarree()
        cmap = sns.light_palette("navy", 12, as_cmap=True)

        fig, ax = plt.subplots(1, 1, figsize=(10, 5),
                               subplot_kw=dict(projection=proj))
        cax = fig.add_axes([0, 0, 0.1, 0.1])
        fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)

        hb = ax.hexbin(lon, lat, gridsize=(50, 15), bins='log', 
                       transform=proj, cmap=cmap)
        
        # This block of code helps correctly size and place the colorbar
        def resize_colorbar(event):
            plt.draw()
            posn = ax.get_position()
            cax.set_position([posn.x0 + posn.width + 0.01, posn.y0, 
                              0.04, posn.height])
        fig.canvas.mpl_connect('resize_event', resize_colorbar)

        ax.coastlines()
        ax.set_global()   

        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=proj)
        ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=proj)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

        plt.colorbar(hb, cax, ax)

        plt.show()

    def n_mu(self):
        import matplotlib.pyplot as plt
        import seaborn as sns

        modes = ['AIT', 'ACC', "MOS", "MBS"]
        df_aero = self._df[['mu%s' % m for m in modes] + ['n%s' % m for m in modes]] 
        df_aero = np.log10(df_aero)

        sel_sh = self._df['lat'] < 0.
        df_aero['hemisphere'] = 'NH'
        df_aero.loc[sel_sh, 'hemisphere'] = 'SH'

        with sns.axes_style("ticks"):
            for mode in modes: 
                g = sns.FacetGrid(df_aero, col="hemisphere")
                g.map(plt.hexbin, 'n'+mode, 'mu'+mode, gridsize=(20, 20), edgecolor='none')

        plt.show()


if __name__ == "__main__":

    qp = QuickPlot()
