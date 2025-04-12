#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Do convex inversion.
"""
import os 
import time
import multiprocessing
import itertools
from argparse import ArgumentParser as ap
import numpy as np

from smfl import do_conv, do_conv_final, golden_spiral_G10


if __name__ == "__main__":
    parser = ap(description="Do convexinversion.")
    parser.add_argument(
        "--lam", type=float, default=None, 
        help="ecliptic longitude")
    parser.add_argument(
        "--beta", type=float, default=None, 
        help="ecliptic latitude")
    parser.add_argument(
        "--Nlam", type=int, default=1, 
        help="number of ecliptic longitude")
    parser.add_argument(
        "--Nbeta", type=int, default=1, 
        help="number of ecliptic latitude")
    parser.add_argument(
        "--N_golden", type=int, default=None, 
        help="Number of poles")
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
        "--resdir", type=str, default="convex_result",
        help="Directory for results of convexinv")
    parser.add_argument(
        "--findir", type=str, default="convex_final",
        help="Directory for results of convexinv")
    parser.add_argument(
        "--final", action="store_true", default=False, 
        help="Create final model with fixed period and pole")
    args = parser.parse_args()
 
    resdir = args.resdir
    findir = args.findir
    os.makedirs(resdir, exist_ok=True)
    os.makedirs(findir, exist_ok=True)
    inpdir = args.inpdir

    lcdir = args.lcdir
    lc = args.lc

 
    # Fixed period and pole
    if args.final:
        do_conv_final(args.lam, args.beta, lc, lcdir, inpdir=inpdir, outdir=findir)

     
    # Starting time 
    t0 = time.time()
    
    if args.N_golden:
        lam, beta = golden_spiral_G10(args.N_golden)
        # Create combinations
        params = [(l, b, lc, lcdir, inpdir, resdir) for l, b in zip(lam, beta)]
        print(f"  Number of iterations N_iter = {len(params)}")
    
    else:
        lam = np.linspace(0, 360, args.Nlam)
        beta = np.linspace(-90, 90, args.Nbeta)
  
        # TODO: accept multiple input
        # Use all *lcs in lcdir
        #lc = f"{lcdir}/*lcs"
        #assert False, "Check the code for multiple lcs."

        # Create combinations
        params = [(l, b, lc, lcdir, inpdir, resdir) 
            for l, b in itertools.product(lam, beta)]
        print(f"  Number of iterations N_iter = {len(params)}")

    n_p = args.N_p
    with multiprocessing.Pool(processes=n_p) as pool:
        pool.starmap(do_conv, params)

    # Ending time 
    t1 = time.time()
    print(f"  Elapsed time : {t1-t0:.1f} s (N_iter = {len(params)})")
