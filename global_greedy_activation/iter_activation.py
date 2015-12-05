#!/usr/bin/env python
""" Compute additional activation details given the results of the global
greedy claculations.

"""
from __future__ import print_function

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import glob
import os
import parcel_model as pm
import pandas as pd

import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger()

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

N_ITER = 4

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("path_to_output", type=str,
                    help="Path to directory containing output CSV files from"
                         " calc_global.py run.")
parser.add_argument("-o", "--output", type=str, default="iter_nact.csv",
                    help="Optional name for output results dataset")
parser.add_argument("-d", "--debug", action="store_true",
                    help="Log extra debugging information")
parser.add_argument("-n", "--n", type=int, default=0,
                    help="Cap for number of values to compute.")

def calc_activation(row, smax):

    V, T, P,                            \
    muAIT, muACC, muMOS, muMBS,         \
    nAIT, nACC, nMOS, nMBS,             \
    nDST01, nDST02, nDST03, nDST04,     \
    nSSLT01, nSSLT02, nSSLT03, nSSLT04, \
    kappaMOS, time, lev, lat, lon       = row

    logger.debug( "      WSUB: %f" % V)
    logger.debug( "       ACC:" + (2*" {}").format(muACC, nACC ))
    logger.debug( "       AIT:" + (2*" {}").format(muAIT, nAIT ))
    logger.debug( "       MBS:" + (2*" {}").format(muMBS, nMBS))
    logger.debug( "       MOS:" + (3*" {}").format(muMOS, nMOS, kappaMOS))
    logger.debug( "      nDST:" + (4*" {}").format(nDST01, nDST02, nDST03, nDST04))
    logger.debug( "     nSSLT:" + (4*" {}").format(nSSLT01, nSSLT02, nSSLT03, nSSLT04))
    logger.info( "    (time, lev, lat, lon) = (%d, %d, %d, %d)" % (time, lev, lat, lon))

    # Setup the aerosol
    aerosol_modes = []
    aer_dict = {}

    aer_dict["ACC"] = pm.AerosolSpecies("ACC",
                            pm.Lognorm(mu=muACC,
                                       sigma=mode_dict['ACC']['sigma'],
                                       N=nACC),
                            kappa=mode_dict['ACC']['kappa'], bins=100)
    aerosol_modes.append(aer_dict["ACC"])

    aer_dict["AIT"] = pm.AerosolSpecies("AIT",
                            pm.Lognorm(mu=muAIT,
                                       sigma=mode_dict['AIT']['sigma'],
                                       N=nAIT),
                            kappa=mode_dict['AIT']['kappa'], bins=100)
    aerosol_modes.append(aer_dict["AIT"])

    aer_dict["MBS"] = pm.AerosolSpecies("MBS",
                            pm.Lognorm(mu=muMBS,
                                       sigma=mode_dict['MBS']['sigma'],
                                       N=nMBS),
                            kappa=mode_dict['MBS']['kappa'], bins=100)
    aerosol_modes.append(aer_dict["MBS"])

    aer_dict["MOS"] = pm.AerosolSpecies("MOS",
                            pm.Lognorm(mu=muMOS,
                                       sigma=mode_dict['MOS']['sigma'],
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

    # Compute number concentration activated
    N_act_tot, N_tot = 0., 0.
    for aer in aerosol_modes:
        N_act, act_frac = pm.lognormal_activation(
            smax, aer.distribution.mu*1e-6, aer.distribution.sigma,
            aer.distribution.N, aer.kappa, T=T
        )
        logger.debug("{} {}".format(aer.species, act_frac))
        N_act_tot += N_act
        N_tot += aer.distribution.N

    logger.debug("----"*20)
    logger.debug("    N_tot = {:f}".format(N_tot))
    logger.debug("N_act_tot = {:f}".format(N_act_tot))
    logger.debug("===="*20)

    return N_act_tot, N_tot

if __name__ == "__main__":

    args = parser.parse_args()

    # Set debug logging if requested
    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug logging enabled")

    logger.debug(args)

    # Check if files in output directory
    logger.info("Reading files in " + args.path_to_output)
    files = sorted(glob.glob(os.path.join(args.path_to_output, "*.csv")))
    if not files:
        logger.critical("WARNING: there were no files found in %s" %
                        args.path_to_output)
    dfs = []
    for f in files:
        logger.info("   " + f)
        dfs.append(pd.read_csv(f, index_col=0))
    df = pd.concat(dfs, ignore_index=True)

    if args.n > 0:
        df = df.ix[:args.n]

    # Function to drive activation calculation
    def map_fn(row):
        results = []
        for i in range(1, N_ITER+1):
            logging.debug("iter %d" % i)
            results.extend(calc_activation(row.ix[:24], row['SMAX_%02d' % i]))
        return results
    results = df.apply(map_fn, axis=1)

    # Create new dataframe
    data = {}
    logging.info("Processing results...")
    for i in range(1, N_ITER+1):
        logging.info("   iter %d" % i)
        data["N_ACT_%02d" % i] = results.apply(lambda s: s[2*(i-1)])
    data["N_TOT"] = results.apply(lambda s: s[1])
    data = pd.DataFrame(data)

    data.to_csv(args.output)



