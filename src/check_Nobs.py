#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check number of observations.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  


if __name__ == "__main__":
    parser = ap(description="Check number of observations.")
    parser.add_argument(
        "lc", help="Preprocessed lightcurves.")
    args = parser.parse_args()
    
    Nobs = 0
    with open(args.lc, "r") as f:
        lines = f.readlines()
        for idx_line, line in enumerate(lines):
            sl = line.split()
            if idx_line == 0:
                Nlc = int(sl[0])
            if len(sl) > 4:
                Nobs += 1

    print(f"N_lc  = {Nlc}")
    print(f"N_obs = {Nobs}")
