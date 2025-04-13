#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Format lightcurves in ALCDEF.
Input lightcurves have time in jd, flux, and fluxerr (for Monte Carlo trials).
Arbitrary keywords can be used with key_jd, key_flux, and key_fluxerr options.
"""
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap

from smfl import format4inv, Ariadnetestdata, save4inv, calc_JPLephem, tbinning


if __name__ == "__main__":
    parser = ap(description="Fotmat lightcurves from ALCDEF")
    parser.add_argument(
        "lc", type=str,
        help="Input from ALCDEF")
    parser.add_argument(
        "out", type=str,
        help="output filename")
    args = parser.parse_args()
   
    jd_list, mag_list, magerr_list = [], [], []
    with open(args.lc, "r") as f:
        lines = f.readlines()

        for l in lines:
            if l[0:5]=="DATA=":
                jd, mag, magerr = l[5:].split("|")
                jd_list.append(float(jd))
                mag_list.append(float(mag))
                magerr_list.append(float(magerr))
            else:
                pass

    # Convert mag to flux
    #   m = -2.5log10(f) + C1
    #   f = 10**(-(m-C1)/2.5)
    #     = C2*10**(-m/2.5)
    flux_list = [10**(-m/2.5) for m in mag_list]
    # Normalize by mean
    C = np.mean(flux_list)
    flux_list = [x/C for x in flux_list]

    # Assume S/N >> 1 
    fluxerr_list = [1-10**(-merr/2.5) for merr in magerr_list]

    
    # TODO: lt correction?
    df = pd.DataFrame(
        dict(jd=jd_list, flux=flux_list, fluxerr=fluxerr_list, mag=mag_list, magerr=magerr_list))

    df.to_csv(args.out, sep=" ")
