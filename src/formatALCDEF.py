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
from astroquery.jplhorizons import Horizons
from astropy import units as u

from smfl import calc_JPLephem


def calc_ltday(obj, code, t_jd):
    """Calculate light traveling time in day.

    Parameters
    ----------
    obj : str
        object name
    code : str
        MPC code
    t_jd : str
        time in jd

    Return
    ------
    lt_day : float
        light traveling time in day
    """
    ast = Horizons(id=obj, location=code, epochs=t_jd)
    eph = ast.ephemerides()
    delta = eph["delta"][0]
    # Calculation is consistent with https://alcdef.org/docs/ALCDEF_Standard.pdf
    lt_day = delta*au2km/c/sec2day
    return lt_day


if __name__ == "__main__":
    parser = ap(description="Fotmat lightcurves from ALCDEF")
    parser.add_argument(
        "lc", type=str,
        help="Input from ALCDEF")
    parser.add_argument(
        "out", type=str,
        help="output filename wo/ extension")
    parser.add_argument(
        "--f_max", type=float, default=False,
        help="Maximum flux to be saved")
    parser.add_argument(
        "--snr_min", type=float, default=False,
        help="Minimum SNR to be saved")
    args = parser.parse_args()
   

    # For light-time correction
    # Light speed 3e5 km/s
    c = 3e5
    # Astronomical unit
    au = 1*u.au
    # in km
    au2km = au.to("km").value
    sec2day = 24.*3600.
    # If lt_day (light traveling time in day) is smaller than this threhold,
    # the mean lt_day is used to perform light traveling time correction.
    lt_th_s = 0.2


    jd_list, mag_list, magerr_list = [], [], []
    nobs_list = []
    # Session date 
    sdate_list = []
    with open(args.lc, "r") as f:
        lines = f.readlines()
        n_obs = 0
        sdate0 = "999"

        for l in lines:
            if l[0:11] == "SESSIONDATE":
                l_sp = l.split("=")
                sdate = l_sp[1].replace("\n", "")
                print(f"SESSION DATE: {sdate} (n_obs={n_obs})")
                if sdate == sdate0:
                    pass
                else:
                    # Update date (n_obs)
                    n_obs += 1
                    sdate_list.append(sdate)
                    sdate0 = sdate
            elif l[0:5]=="DATA=":
                jd, mag, magerr = l[5:].split("|")
                jd_list.append(float(jd))
                mag_list.append(float(mag))
                magerr_list.append(float(magerr))
                nobs_list.append(n_obs)
            elif l[0:7] == "ENDDATA":
                # Update n_obs ? -> No. It is better to use SESSIONDATE
                pass
            elif l[0:7] == "MPCCODE":
                l_sp = l.split("=")
                code = l_sp[1].replace("\n", "")
                print(f"MPCCODE: {code} (n_obs={n_obs})")
            else:
                pass

    df = pd.DataFrame({"jd": jd_list, "mag":mag_list, "magerr":magerr_list, "n_obs": nobs_list})

    df_list = []
    key_mag, key_magerr = "mag", "magerr"
    df["lt_day"] = 0
    for nobs in sorted(list(set(nobs_list))):
        df_temp = df[df["n_obs"] == nobs]

        # Set mean mag to 0
        df_temp[key_mag] = df_temp[key_mag] - np.mean(df_temp[key_mag])
        # F    = 10**(-0.4*mag)
        # Ferr = 0.4*log_e(10)*F*magerr
        df_temp["flux_cor"] = 10**(-0.4*df_temp[key_mag])
        df_temp["fluxerr_cor"] = 0.4*np.log(10)*df_temp["flux_cor"]*df_temp[key_magerr]

        # Lighttime correction
        t_jd0 = np.min(df_temp["jd"])
        t_jd1 = np.max(df_temp["jd"])
        lt_day0 = calc_ltday("2015 BY310", code, t_jd0)
        lt_day1 = calc_ltday("2015 BY310", code, t_jd1)
        print(f"lt_day0: {lt_day0:.2f} d = {lt_day0*24.:.2f} hr = {lt_day0*24.*3600.:.2f} s")
        print(f"lt_day1: {lt_day1:.2f} d = {lt_day1*24.:.2f} hr = {lt_day1*24.*3600.:.2f} s")
        lt_s = (lt_day1 - lt_day0)*24.*3600.
        assert lt_s < lt_th_s
        
        t_jd_mean = (t_jd0 + t_jd1)/2.
        lt_day_mean = calc_ltday("2015 BY310", code, t_jd_mean)
        df_temp["lt_day"] = lt_day_mean
        df_list.append(df_temp)

    df = pd.concat(df_list)

    # Do light-time correction
    df["t_jd_ltcor"] = df["jd"] - df["lt_day"]
    
    out = f"{args.out}.txt"
    df.to_csv(out, sep=" ", index=False)

    # Save separately
    for nobs, fname in zip(sorted(list(set(nobs_list))), sdate_list):
        df_temp = df[df["n_obs"] == nobs]
        # Remove outliers
        if args.f_max:
            df_temp = df_temp[df_temp["flux_cor"] < args.f_max]
        if args.snr_min:
            df_temp = df_temp[df_temp["flux_cor"]/df_temp["fluxerr_cor"] > args.snr_min]
        out = f"{args.out}_{fname}.txt"
        df_temp.to_csv(out, sep=" ", index=False)
