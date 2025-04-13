#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Format lightcurves downloaded from ALCDEF (https://alcdef.org/).

Output file has such as "jd (time in jd)", "flux", and "fluxerr" (for Monte Carlo trials).
Arbitrary keywords can be used with key_jd, key_flux, and key_fluxerr options.
After using this script, you can format lightcurves for convex inversion with "format4convexinv.py".
"""
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap

from smfl import calc_JPLephem


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

    df = pd.DataFrame({"jd": jd_list, "mag":mag_list, "magerr":magerr_list})

    # Set mean mag to 0
    key_mag, key_magerr = "mag", "magerr"
    df[key_mag] = df[key_mag] - np.mean(df[key_mag])
    # F    = 10**(-0.4*mag)
    # Ferr = 0.4*log_e(10)*F*magerr
    df["flux_cor"] = 10**(-0.4*df[key_mag])
    df["fluxerr_cor"] = 0.4*np.log(10)*df["flux_cor"]*df[key_magerr]

    ## TODO: lt correction?
    #df = pd.DataFrame(
    #    dict(jd=jd_list, flux=flux_list, fluxerr=fluxerr_list, mag=mag_list, magerr=magerr_list))

    df.to_csv(args.out, sep=" ", index=False)
