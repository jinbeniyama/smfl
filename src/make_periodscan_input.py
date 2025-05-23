#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make input files of period_scan.
p = 0.8 (for coarse search) and p=0.1 (fine search) are good choices.

Note from convexinv_doc.pdf
---------------------------
The first line of the input period scan gives the initial period, the final period, and the
coeffient p of the period step. The interval is scanned with the step p ∆P (always set p < 1,
a recommended value is p = 0.8). Other input parameters have the same meaning as in the
input par file. The last line gives the minimum number of iterations if the iteration stop
condition is smaller than one.


- CRW (Convexity regularization weight), 0.1 by default:
  This is often needed to keep the shape formally convex.
  A typical value is 0.1, but it may often have to be orders of magnitude larger or smaller.
  There is a dark facet area that makes the whole set of facets convex. 
  Always try to put it below 1% of the total area by increasing the convexity 
  regularization parameter. 
  If you cannot get a good fit with a small dark area and 
  if the fit significantly improves when the dark facet is large, 
  it means that there is likely to be an albedo variation over the surface.

- ISC (Iteration stop condition), 50 by default:
  If it is an integral number higher than one, 
  then it is the number of iteration steps in the Levenberg-Marquardt loop. 
  If it is lower than one, then it is the smallest difference in rms 
  deviation between two subsequent steps – when the steps have smaller 
  difference, the iteration loop is stopped.

- MNI (Minimum number of iterations), 10 by default:

Template
--------
0.0925 0.09277 0.1	period start - end - interval coeff.
0.1			convexity weight
6 6			degree and order of spherical harmonics
8			no. of rows
0.5	0		scattering parameters
0.1	0
-0.5	0
0.1	0
50	  		iteration stop condidion
10			minimum number of iterations (only if the above value < 1)
--------
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == "__main__":
    parser = ap(description="Create input of convexinv with arbitary settings")
    parser.add_argument(
        "--P0_hr", type=float, default=0, 
        help="Minimum period in hour")
    parser.add_argument(
        "--P1_hr", type=float, default=0, 
        help="Minimum period in hour")
    parser.add_argument(
        "--dt_hr", type=float, default=24., 
        help="Observation arc in hour")
    parser.add_argument(
        "--interval", type=float, default=0.8, 
        help="Interval coefficient")
    parser.add_argument(
        "--deg_harmonics", type=int, default=6, 
        help="Degree of spherical harmonics")
    parser.add_argument(
        "--ord_harmonics", type=int, default=6, 
        help="Order of spherical harmonics")
    parser.add_argument(
        "--CRW", type=float, default=0.1, 
        help="Convexity regularization weight")
    parser.add_argument(
        "--ISC", type=float, default=50, 
        help="Iteration stop condition")
    parser.add_argument(
        "--MNI", type=int, default=10, 
        help="Minimum number of iterations")
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
        "--out", type=str, default="input_ps",
        help="Output file name")
    args = parser.parse_args()
   
    P0 = args.P0_hr
    P1 = args.P1_hr
    eps = args.interval
    N = args.deg_harmonics
    M = args.ord_harmonics
    CRW = args.CRW
    ISC = args.ISC
    MNI = args.MNI

    N_per = 2*args.dt_hr*(P1-P0)/P0/P1/eps
    print(f"Number of periods: N_per = {N_per:.1f}")
    
    with open(args.out, "w") as f:
        f.write(f"{P0} {P1} {eps} \n")
        # Convex weight
        f.write(f"{CRW} \n")
        # Degree and order of spherical harmonics
        f.write(f"{N} {M}\n")
        f.write("8\n")
        f.write(f"{args.a} 0\n")
        f.write(f"{args.d} 0\n")
        f.write(f"{args.k} 0\n")
        f.write(f"{args.c} 0\n")
        f.write(f"{ISC}\n")
        f.write(f"{MNI}")
        
