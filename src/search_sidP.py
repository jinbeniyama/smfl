#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Search siderial period (sidP) with period_scan.
"""
from argparse import ArgumentParser as ap
import numpy as np
import subprocess


if __name__ == "__main__":
    parser = ap(description="Search sidereal period with period_scan")
    parser.add_argument(
        "lc", type=str, nargs="*",
        help="lightcurves (multiple lcs are acceptable only for MC technique)")
    parser.add_argument(
        "--inp", type=str, default="input_ps", 
        help="input file of period_scan")
    parser.add_argument(
        "--out", type=str, 
        help="output file of period_scan")
    args = parser.parse_args()


    if args.out:
        out = args.out
    else:
        # Save in current directory
        out = f"{args.inp.split('/')[-1]}_out"

    for idx,lc in enumerate(args.lc):
        cmd = f"cat {lc} | period_scan -v {args.inp} {out}_{idx+1:04d}"
        print(f"{idx+1:04d} : {cmd}")
        subprocess.run(cmd, shell=True)
