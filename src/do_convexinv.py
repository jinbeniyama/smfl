#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Do convex inversion.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import subprocess


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
   
    # Initial pole direction
    lam0  = 0
    beta0 = 0
  
    lam = np.linspace(0, 360, args.Nlam)
    beta = np.linspace(-90, 90, args.Nbeta)
  
    if args.lcdir:
        # Use all *lcs in lcdir
        lcdir = args.lcdir
        lc = f"{lcdir}/*lcs"
    else:
        lc = args.lc
  
    for l in lam:
        for b in beta:
            inp = f"input_ci_{int(l)}_{int(b)}"
            outarea = f"outarea_ci_{int(l)}_{int(b)}"
            outpar = f"outpar_ci_{int(l)}_{int(b)}"
            outlc = f"outlcs_ci_{int(l)}_{int(b)}"
            # Include chi2
            res = f"res_ci_{int(l)}_{int(b)}"
            cmd = (
                f"cat {lc} | convexinv -v -o {outarea} "
                f"-p {outpar} {inp} {outlc} > {res}"
                )
            subprocess.run(cmd, shell=True)
