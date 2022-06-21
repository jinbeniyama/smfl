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
"""
from argparse import ArgumentParser as ap
import numpy as np

from smlf import format4inv, save4inv


if __name__ == "__main__":
  parser = ap(description="Fotmat lightcurves for convex inversion.")
  parser.add_argument(
    "csv", nargs="*",
    help="Results of light curve observation")
  parser.add_argument(
    "--jpl", nargs="*", required=True,
    help="JPL/HORIZONS result")
  parser.add_argument(
    "--absflux", action="store_true", default=False,
    help="Use absolute value of flux")
  parser.add_argument(
    "--N_mc", type=int, default=1, 
    help="number of trials using Monte Carlo technique")
  parser.add_argument(
    "--out", default="out_lcs", 
    help="output filename")
  args = parser.parse_args()
 
  assert len(args.csv)==len(args.jpl), "phot/jpl res should be the same dim"
  
  result = []
  # Extract object info
  for csv,jpl in zip(args.csv, args.jpl):
      df_phot = format4inv(csv, jpl, args.absflux)
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
      save4inv(result, args.absflux, random, out)
