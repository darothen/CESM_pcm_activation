
from __future__ import print_function

from datetime import datetime
from numpy import pi, exp, log, timedelta64, max, min
import pandas as pd
import xray

import os
import sys

from argparse import ArgumentParser
parser = ArgumentParser(
    description="Compute aerosol size distribution properties.")
parser.add_argument
parser.add_argument("in_file",
                    help="MARC aerosol output file to process")
parser.add_argument("-o", "--output",
                    help="Name of output file")
parser.add_argument("-y", "--years-offset", type=int, default=0,
                    help="Offset for aligning input file times")
parser.add_argument("-a", "--aer-props", default="aer_props.csv",
                    help="Path to aerosol properties table")
parser.add_argument("--parallel", action="store_true",
                    help="Pre-chunk the data so that xray uses dask")


def _timestamp():
    """ Timestamp the current date/time """
    now = datetime.now()
    return datetime.strftime(now, "%a %b %-d %H:%M:%S %Y")

if __name__ == "__main__":

    args = parser.parse_args()
    print(" ".join(sys.argv))

    # File I/O setup
    in_fn_base = os.path.basename(args.in_file)
    if args.output is None:
        base, ext = os.path.splitext(in_fn_base)
        out_fn = base + ".aerosol_dists" + ext
    else:
        out_fn = args.output

    # Open and read file
    print("Opening %s" % args.in_file)
    data = xray.open_dataset(args.in_file, decode_times=False)
    if args.parallel:
        data = data.chunk({'time': 10})

    # Fix time alignment
    data_coords = data.coords.to_dataset()
    data_coords.time.values += args.years_offset*365.
    data = xray.decode_cf(data)
    itimes = data.time.to_index()
    itimes = itimes + timedelta64(9, 'Y')
    data = data.assign_coords(time=itimes)

    # Select timeslice 3 - end
    data = data.isel(time=slice(58, None))

    Rd = 287. #: Gas constant, dry air  [J/kg]

    # Compute thermodynamic state variables at each grid cell
    P = data.hyam*data.P0 + data.hybm*data.PS
    w = data.Q/(1.0 - data.Q)
    Tv = data['T']*(1.0 + 0.61*w)

    rho = P/(Rd*Tv)
    data['rho'] = rho
    data['rho'].attrs['units'] = "kg/m3"
    data['rho'].attrs['long_name'] = 'air density'


    _TAB = "  "

    MARC_MODES = ["OC", "BC", "MOS", "MBS", "NUC", "AIT", "ACC"]
    MARC_NAMES = ["organic carbon", "black carbon", "mixed OC-sulfate",
                  "mixed BC-sulfate", "nucleation sulfate", "aitken sulfate",
                  "accumulation sulfate"]

    SSLT_MODES = ["SSLT%02d" % i for i in range(1, 5)]
    SSLT_NAMES = ["sea salt %d" % i for i in range(1, 5)]

    DST_MODES = ["DST%02d" % i for i in range(1, 5)]
    DST_NAMES = ["dust %d" % i for i in range(1, 5)]

    aer_props = pd.read_csv("aer_props.csv", index_col=0)

    # Number concentration
    print("Number concentration...")
    for mode in MARC_MODES:
        print(_TAB + mode)
        data["n"+mode] = data["n"+mode]*rho*1e-6
        data["n"+mode].attrs['units'] = "1/cm3"

    for mode in SSLT_MODES + DST_MODES:
        print(_TAB + mode)
        num_to_mass = aer_props.ix[mode].num_to_mass
        data["n"+mode] = data[mode]*rho*num_to_mass*1e-6
        data["n"+mode].attrs['units'] = "1/cm3"

    # Mass concentrations, kg/kg -> ng/m3
    print("Mass concentration...")
    for mode in MARC_MODES + ['OIM', 'BIM']:
        print(_TAB + mode)
        data["m"+mode] = data["m"+mode]*rho*1e9
        data["m"+mode].attrs['units'] = "ng/cm3"

    for mode in SSLT_MODES + DST_MODES:
        print(_TAB + mode)
        data["m"+mode] = data[mode]*rho*1e9
        data["m"+mode].attrs['units'] = "1/cm3"

        # Can now safely delete old single-moment modes
        del data[mode]

    # Mixed-mode hygroscopicities
    print("Mixed-mode hygroscopicities")

    def density(mode):
        if mode == 'MOS':
            return data['rhoMOS']
        else:
            return aer_props.ix[mode].density

    def kappa(mode):
        if mode == 'MOS':
            return data['kappaMOS']
        else:
            return aer_props.ix[mode].hygro

    data['rhoMOS'] = ( data.mOIM*density('OC') +
                       (data.mMOS - data.mOIM)*density('ACC') ) / data.mMOS
    # Limit min/max
    min_density = min([density('OC'), density('ACC')])
    max_density = max([density('OC'), density('ACC')])
    data['rhoMOS'].values[data.rhoMOS.values > max_density] = max_density
    data['rhoMOS'].values[data.rhoMOS.values < min_density] = min_density
    data['rhoMOS'].attrs['units'] = 'kg/m3'
    data['rhoMOS'].attrs['long_name'] = \
        'avg. particle density in mixed organic/sulfate mode'

    data['kappaMOS'] = ( kappa('OC')*data.mOIM*1e-9/density('OC') +
                         kappa('ACC')*(data.mMOS - data.mOIM)*1e-9/density('ACC')
                        ) / (data.mMOS*1e-9 / data.rhoMOS)
    min_kappa = min([kappa('OC'), kappa('ACC')])
    max_kappa = max([kappa('OC'), kappa('ACC')])
    data['kappaMOS'].values[data.kappaMOS.values > max_kappa] = max_kappa
    data['kappaMOS'].values[data.kappaMOS.values < min_kappa] = min_kappa
    data['kappaMOS'].attrs['units'] = ""
    data['kappaMOS'].attrs['long_name'] = \
        "hygroscopicity of mixed organic/sulfate mode"

    # Modal mean sizes
    def sigma(mode):
        return aer_props.ix[mode].dispersion

    def calc_mu(mode):
        num = data['n' + mode]
        mass = data['m' + mode]

        density = aer_props.ix[mode].density
        sigma = aer_props.ix[mode].dispersion

        norm = ((mass*1e-9)*(3./4./pi)/density/num)
        exp_fact = exp(-4.5*(log(sigma)**2))
        mu_cubed = norm*exp_fact
        mu = mu_cubed**(1./3.)
        mu *= 1e4 # -> convert to micron

        return mu

    print("Modal mean sizes")
    for mode, name in zip(MARC_MODES, MARC_NAMES):
        print(_TAB + mode)
        data['mu' + mode] = calc_mu(mode)
        data['mu' + mode].attrs['units'] = "micron"
        data['mu' + mode].attrs['long_name'] = \
            "geo. mean radius, %s" % name

    # Append script info to history
        his_kws = dict(time=_timestamp(), command=" ".join(sys.argv))

    if "history" in data.attrs:
        history = data.attrs['history']
        history = "\n{time:s}:  {command:s}".format(**his_kws) + history
        data.attrs['history'] = history
    else:
        data.attrs['history'] = "{time:s}:  {command:s}".format(**his_kws)

    print("Writing output (could take a while)...")
    data.to_netcdf(out_fn)
