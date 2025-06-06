#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make input files of convexinv in outdir/.

Note from convexinv_doc.pdf
---------------------------
Start with, e.g., a = 0.5, d = 0.1, k = −0.5. 
If you use only relative lightcurves, no solar phase function is needed 
(i.e., relative fit with almost constant solar phase for each lightcurve), 
then use fixed values of phase function parameters with zero amplitude 
and slope. 
Fitting the phase parameters using relative data leads to divergence. 
The parameter c has usually only little effect on the solution 
so you can fix it at, e.g., c = 0.1. 

Template
--------
100		0	inital lambda [deg] (0/1 - fixed/free)
100		0	initial beta [deg] (0/1 - fixed/free)
0.0143		1	inital period [hours] (0/1 - fixed/free)
0			zero time [JD]
0			initial rotation angle [deg]
0.1			convexity regularization
6 6			degree and order of spherical harmonics expansion
8			number of rows
0.5		0	phase funct. param. 'a' (0/1 - fixed/free)
0.1		0	phase funct. param. 'd' (0/1 - fixed/free)
-0.5		0	phase funct. param. 'k' (0/1 - fixed/free)
0.1		0	Lambert coefficient 'c' (0/1 - fixed/free)
50			iteration stop condition
--------

"""
import os 
from argparse import ArgumentParser as ap
import numpy as np

from smfl import golden_spiral_G10, mean_angular_spacing

if __name__ == "__main__":
    parser = ap(description="Create input of convexinv with arbitary settings")
    parser.add_argument(
        "--sidP", type=float, default=0, 
        help="sidereal period in hour")
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
        "-a", type=float, default=0.5, 
        help="phase funct. param a (amplitude)")
    parser.add_argument(
        "-d", type=float, default=0.1, 
        help="phase funct. param d (width)")
    parser.add_argument(
        "-k", type=float, default=-0.5, 
        help="phase funct. param k (slope)")
    parser.add_argument(
        "-c", type=float, default=0.1, 
        help="Lambert coefficient c")
    parser.add_argument(
        "--CRW", type=float, default=0.1, 
        help="Convexity regularization weight")
    parser.add_argument(
        "--ISC", type=float, default=50, 
        help="Iteration stop condition")
    parser.add_argument(
        "--fixP", action="store_true", default=False, 
        help="Fix sidereal period")
    parser.add_argument(
        "--fixpole", action="store_true", default=False, 
        help="Fix orientations")
    parser.add_argument(
        "--inpdir", type=str, default="convex_input",
        help="Directory for input file of convexinv")
    parser.add_argument(
        "--out", type=str, default=False,
        help="Output file name")
    args = parser.parse_args()
   
    outdir = args.inpdir
    os.makedirs(outdir, exist_ok=True)

    CRW = args.CRW
    ISC = args.ISC
    
    if args.lam:
        lam = [args.lam]
        beta = [args.beta]
    elif args.N_golden:
        lam, beta = golden_spiral_G10(args.N_golden)
    else:
        lam = np.linspace(0, 360, args.Nlam)
        beta = np.linspace(-90, 90, args.Nbeta)
    
    # siderial period in hour
    if args.sidP == 0:
        sidP = 0
    else:
        sidP = args.sidP
  
    # Whether fix sidP
    if args.fixP:
        fix_P = 0
    else:
        fix_P = 1
  
    # Whether fix pole orientaions
    if args.fixpole:
        fix_beta   = 0
        fix_lambda = 0
    else:
        fix_beta   = 1
        fix_lambda = 1
    
    if args.N_golden:
        for l, b in zip(lam, beta):
            # This is for N(lam)==N(beta)==1
            if args.out:
                out = args.out
            else:
                out = f"input_ci_{int(l)}_{int(b)}"
            out = os.path.join(outdir, out)
            with open(out, "w") as f:
                f.write(f"{l} {fix_lambda}\n")
                f.write(f"{b} {fix_beta}\n")
                f.write(f"{sidP} {fix_P}\n")
                # Zero time
                f.write(f"0\n")
                # Initial rotation angle
                f.write(f"0\n")
                # convexity regularization      
                f.write(f"{CRW}\n")
                # degree and order of spherical harmonics expansion
                f.write(f"6 6\n")
                # number of rows        
                f.write(f"8\n")
                # phase funct. param. 'a' (0/1 - fixed/free)  
                f.write(f"{args.a} 0\n")
                # phase funct. param. 'd' (0/1 - fixed/free) 
                f.write(f"{args.d} 0\n")
                # phase funct. param. 'k' (0/1 - fixed/free) 
                f.write(f"{args.k} 0\n")
                # Lambert coefficient 'c' (0/1 - fixed/free)      
                f.write(f"{args.c} 0\n")
                # iteration stop condition
                f.write(f"{ISC}")
    else:
        for l in lam:
            for b in beta:
                # This is for N(lam)==N(beta)==1
                if args.out:
                    out = args.out
                else:
                    out = f"input_ci_{int(l)}_{int(b)}"
                out = os.path.join(outdir, out)
                with open(out, "w") as f:
                    f.write(f"{l} {fix_lambda}\n")
                    f.write(f"{b} {fix_beta}\n")
                    f.write(f"{sidP} {fix_P}\n")
                    # Zero time
                    f.write(f"0\n")
                    # Initial rotation angle
                    f.write(f"0\n")
                    # convexity regularization      
                    f.write(f"{CRW}\n")
                    # degree and order of spherical harmonics expansion
                    f.write(f"6 6\n")
                    # number of rows        
                    f.write(f"8\n")
                    # phase funct. param. 'a' (0/1 - fixed/free)  
                    f.write(f"{args.a} 0\n")
                    # phase funct. param. 'd' (0/1 - fixed/free) 
                    f.write(f"{args.d} 0\n")
                    # phase funct. param. 'k' (0/1 - fixed/free) 
                    f.write(f"{args.k} 0\n")
                    # Lambert coefficient 'c' (0/1 - fixed/free)      
                    f.write(f"{args.c} 0\n")
                    # iteration stop condition
                    f.write(f"{ISC}")
