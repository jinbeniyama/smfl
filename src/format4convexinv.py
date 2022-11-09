#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Format lightcurves for convexinv.
Calculate x_sun, y_sun ,z_sun (sun cartesian coodinate from the object)
and x_earth, y_earth, z_earth (earth cartesian coordinate from the object).

Input lightcurves have time in jd, flux, and fluxerr (for Monte Carlo trials).
Arbitrary keywords can be used with key_jd, key_flux, and key_fluxerr options.

Before execute this script, query JPL/HORIZONS and obtain quantities below.
  'Date__(UT)__HR:HM', 'R.A._____(ICRF)_____DEC', 'R.A.__(a-apparent)__DEC',
  'hEcl-Lon', 'hEcl-Lat', 'r', 'rdot', 'delta', 'deldot'

Read a documentation in DAMIT website about the format of output.


Test: Conpare the output "ariadne_test", "test_lcs_rel" in DAMIT script,
      and lc.txt in "https://astro.troja.mff.cuni.cz/projects/damit/LightCurves
      /exportAllForAsteroid/122/plaintext/A122.lc.txt". 
      Flux values should be ignored.
      (JD is slightly dirrecnt between test_lcs_rel and lc.txt....)

      format4convexinv.py --test

"""
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap

from smfl import format4inv, Ariadnetestdata, save4inv, calc_JPLephem, tbinning


if __name__ == "__main__":
    parser = ap(description="Fotmat lightcurves for convex inversion.")
    parser.add_argument(
        "--res", nargs="*",
        help="Results of light curve observation")
    parser.add_argument(
        "--jpl", nargs="*",
        help="JPL/HORIZONS result")
    parser.add_argument(
        "--absflux", action="store_true", default=False,
        help="Use absolute value of flux")
    parser.add_argument(
        "--N_mc", type=int, default=1, 
        help="number of trials using Monte Carlo technique")
    parser.add_argument(
        "--key_jd", default="t_jd_ltcor", 
        help="Keyword of time in JD")
    parser.add_argument(
        "--key_flux", default="flux", 
        help="Keyword of flux")
    parser.add_argument(
        "--key_fluxerr", default="fluxerr", 
        help="Keyword of flux uncertainty")
    parser.add_argument(
        "--tbin", type=float, default=None, 
        help="width of time bin")
    parser.add_argument(
        "--unit_t", default="s", 
        help="Unit of time for binning")
    parser.add_argument(
        "--out", default="out_lcs", 
        help="output filename")
    parser.add_argument(
        "--test", action="store_true", default=False, 
        help="Test for Ariadne")
    args = parser.parse_args()
   
  
    # Test
    # Need JPL ephemeris before this script.
    if args.test:
        # Save Ariadnetestdata
        df_phot = Ariadnetestdata()
        csv_Ariadne = "Ariadnetest_data.csv"
        df_phot.to_csv(csv_Ariadne, sep=" ", index=False)
  
        # Save Ariadne JPL ephemeris
        # date0, date1 = "1965-05-01", "1965-05-02"
        # step, code = "1m", "381"
        # eph = calc_JPLephem("Ariadne", date0, date1, step, code)
        # jpl_Ariadne = "Ariadnetestephem.csv"
        # eph.write(jpl_Ariadne, format="pandas.csv", sep=" ")
  
        # Save Ariadne JPL ephemeris by hand
        # Leiden Station, Johannesburg (observatory) [code: 081]3 
        # Reference : vanHouten-Groeneveld et al. 1979
        jpl_Ariadne = "Ariadnetest_jpl.csv"

        df_phot = format4inv(csv_Ariadne, jpl_Ariadne, args.key_jd)
        out = f"test_Ariadne"
        save4inv(
            [df_phot], args.absflux, random, args.key_jd, args.key_flux,
            args.key_fluxerr, out)
        sys.exit()
  
    assert len(args.res)==len(args.jpl), "phot/jpl res should be the same dim"
  
    result = []
    # Extract object info
    for csv,jpl in zip(args.res, args.jpl):
        df_phot = pd.read_csv(csv, sep=" ")
        # Binning with time
        if args.tbin:
            df_phot = tbinning(
                df_phot, args.tbin, key_t=args.key_jd, unit_t="d", 
                key_flux=args.key_flux, key_fluxerr=args.key_fluxerr)
      
        df_phot = format4inv(df_phot, jpl, args.key_jd)
        result.append(df_phot)
  
    # Set seed 
    seed = 0
    np.random.seed(seed)
  
    # If n_mc > 0, save N_mc outputs with Monte Carlo technique
    for n_mc in range(args.N_mc):
        if args.N_mc == 1:
            out = args.out
        else:
            out = f"{args.out}_{n_mc+1:04d}"
        print(f"  N_mc : {n_mc+1}, output : {out}")
  
        # Output results
        if n_mc > 0:
            random = True
        else:
            random = False
        save4inv(
            result, args.absflux, random, args.key_jd, args.key_flux,
            args.key_fluxerr, out)
