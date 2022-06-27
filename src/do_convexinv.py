#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Do convex inversion.
"""
import os 
import time
import multiprocessing
import itertools
from argparse import ArgumentParser as ap
import numpy as np

from smfl import do_conv


if __name__ == "__main__":
    parser = ap(description="Do convexinversion.")
    parser.add_argument(
        "--Nlam", type=int, default=1, 
        help="number of ecliptic longitude")
    parser.add_argument(
        "--Nbeta", type=int, default=1, 
        help="number of ecliptic latitude")
    parser.add_argument(
        "--lc", type=str, default=None, 
        help="preprocessed lightcurve")
    parser.add_argument(
        "--lcdir", type=str, default=None, 
        help="lightcurve directory")
    args = parser.parse_args()
 
    # Starting time 
    t0 = time.time()

    # Initial pole direction
    lam0  = 0
    beta0 = 0
  
    lam = np.linspace(0, 360, args.Nlam)
    beta = np.linspace(-90, 90, args.Nbeta)
  
    if args.lcdir:
        # Use all *lcs in lcdir
        lcdir = args.lcdir
        lc = f"{lcdir}/*lcs"
        assert False, "Check the code for multiple lcs."
    else:
        lc = args.lc

    # Create combinations
    params = [(l, b, lc) for l, b in itertools.product(lam, beta)]
    print(f"  Number of iterations N_iter = {len(params)}")

    n_p = 4
    with multiprocessing.Pool(processes=n_p) as pool:
        pool.starmap(do_conv, params)

    # Ending time 
    t1 = time.time()
    print(f"  Elapsed time : {t1-t0:.1f} s (N_iter = {len(params)})")
