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
        "--N_p", type=int, default=4, 
        help="number of processes")
    parser.add_argument(
        "--lc", type=str, default=None, 
        help="preprocessed lightcurve")
    parser.add_argument(
        "--lcdir", type=str, default=".", 
        help="lightcurve directory")
    parser.add_argument(
        "--inpdir", type=str, default="convex_input",
        help="Directory for input file of convexinv")
    parser.add_argument(
        "--outdir", type=str, default="convex_result",
        help="Directory for input file of convexinv")
    args = parser.parse_args()
 
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    inpdir = args.inpdir
    lcdir = args.lcdir

    # Starting time 
    t0 = time.time()

    lam = np.linspace(0, 360, args.Nlam)
    beta = np.linspace(-90, 90, args.Nbeta)
  
    # TODO: accept multiple input
    # Use all *lcs in lcdir
    #lc = f"{lcdir}/*lcs"
    #assert False, "Check the code for multiple lcs."
    lc = args.lc

    # Create combinations
    params = [(l, b, lc, lcdir, inpdir, outdir) 
        for l, b in itertools.product(lam, beta)]
    print(f"  Number of iterations N_iter = {len(params)}")

    n_p = args.N_p
    with multiprocessing.Pool(processes=n_p) as pool:
        pool.starmap(do_conv, params)

    # Ending time 
    t1 = time.time()
    print(f"  Elapsed time : {t1-t0:.1f} s (N_iter = {len(params)})")
