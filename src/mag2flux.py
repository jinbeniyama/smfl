#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert magnitude to intensity for convexinv.
"""

import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd


if __name__ == "__main__":
    parser = ap(description="Convert magnitude to flux")
    parser.add_argument(
        "res", type=str,
        help="photometric light curve in mag")
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
        "--key_mag", default="mag_red_phase", 
        help="Keyword of mag")
    parser.add_argument(
        "--key_magerr", default="magerr", 
        help="Keyword of magnitude uncertainty")
    parser.add_argument(
        "--out", default="cor_lc.txt", 
        help="output filename")
    args = parser.parse_args()
    
    key_jd = args.key_jd
    key_mag = args.key_mag
    key_magerr = args.key_magerr
    key_flux = args.key_flux
    key_fluxerr = args.key_fluxerr

    df = pd.read_csv(args.res, sep=" ")
    assert set([key_jd, key_mag, key_magerr]) <= set(df.columns.tolist())
    
    # Set mean mag to 0
    df[key_mag] = df[key_mag] - np.mean(df[key_mag])
    # F    = 10**(-0.4*mag)
    # Ferr = 0.4*log_e(10)*F*magerr
    df["flux_cor"] = 10**(-0.4*df[key_mag])
    df["fluxerr_cor"] = 0.4*np.log(10)*df["flux_cor"]*df[key_magerr]

    df.to_csv(args.out, sep=" ")
