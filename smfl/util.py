#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Handle light curves for convexinv.

Query JPL/HORIZONS and obtain quantities below.
  'Date__(UT)__HR:HM', 'R.A._____(ICRF)_____DEC', 'R.A.__(a-apparent)__DEC',
  'hEcl-Lon', 'hEcl-Lat', 'r', 'rdot', 'delta', 'deldot'

Calculate x_sun, y_sun ,z_sun (sun cartesian coodinate from the object)
and x_earth, y_earth, z_earth (earth cartesian coordinate from the object).
"""
import os 
import subprocess
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
from scipy import interpolate
import datetime
from astroquery.jplhorizons import Horizons
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy
from astropy.time import Time


def time_keeper(func):
    """
    Decorator to measure time.
    """
    # To take over docstring etc.
    @wraps(func)
    def wrapper(*args, **kargs):
        t0 = time.time()
        result = func(*args, **kargs)
        t1 = time.time()
        t_elapse = t1 - t0
        print(f"[time keeper] t_elapse = {t_elapse:.03f} s (func :{func.__name__})")
        return result
    return wrapper


def format4inv(lc, jpleph, key_jd):
    """
    Format lightcurves for convex inversion.

    # ToDo: lc and jpleph should be array-like 

    Parameters
    ----------
    lc : str
        filename of lightcurve data
    jpleph : str
        filename of ephemeris from JPL
    key_jd : str
        keyword of time in JD

    Return
    ------
    df_format : pandas.DataFrame
        formatted data
    """
    # Create result = [df_phot1, df_phot2, ...] with jpl info.
    # create jpl dataframe
    utc_list, ra_list, dec_list = [], [], []
    hecl_lon, hecl_lat, r_list, delta_list = [], [], [], []
    jd_list, mjd_list = [], []
    x_sun_list, y_sun_list, z_sun_list = [], [], []
    x_earth_list, y_earth_list, z_earth_list = [], [], []
    with open(jpleph, "r") as f:
        # skip header of jpl ephem 
        n_header = 84
        #n_use = 121
        #f = f.readlines()[n_header:n_header+n_use]
        f = f.readlines()[n_header:]
        for line in f:
          # ok
          utc = line[1:18]
          other = line[18:]
          # Split by space
          other = other.split(" ")
          # Replace "\n" with ""
          other = [x.replace("\n", "") for x in other]
          other = [x for x in other if x != ""]
          if len(other)!=19:
              continue

          utc_list.append(utc)
          ra    = f"{other[1]} {other[2]} {other[3]}"
          dec   = f"{other[4]} {other[5]} {other[6]}"
          beta  = f"{other[13]}"
          lam   = f"{other[14]}"
          r     = f"{other[15]}"
          delta = f"{other[17]}"

          ra_list.append(ra)
          dec_list.append(dec)
          hecl_lon.append(beta)
          hecl_lat.append(lam)
          r_list.append(r)
          delta_list.append(delta)

          # jd
          ## convert utc (2020-Nov-30 16:00) to '1999-01-01T00:00:00.123'
          utc = datetime.datetime.strptime(utc, "%Y-%b-%d %H:%M")
          utc = datetime.datetime.strftime(utc, "%Y-%m-%dT%H:%M:%S")
          t = Time(utc, format='isot', scale='utc')
          jd_list.append(t.jd)
          mjd_list.append(t.mjd)
 
          # the Sun x, y, z
          ## use cartisian method
          c_eclip = SkyCoord(
              beta, lam, r, frame="heliocentricmeanecliptic", 
              unit=(u.deg, u.deg, u.au)
              )
          #print("The Sun ========================================================")
          #print("eclip lon, lat, r from the Sun to the asteroid")
          #print(f"{c_eclip.lon}, {c_eclip.lat}, {c_eclip.distance}")

          c_eclip = c_eclip.cartesian
          ## minus means from object
          #print("converted cartesian")
          x_sun = -c_eclip.x.value
          y_sun = -c_eclip.y.value
          z_sun = -c_eclip.z.value
          #print(f"{x_sun}, {y_sun}, {z_sun}")
          #print("")

          #print("The Earth ======================================================")
          # not icrs (), but gcrs(Geocentric Celestial Reference System)
          c_radec = SkyCoord(
              ra=ra, dec=dec, distance=delta, frame="gcrs", unit=(u.hourangle, u.deg, u.au),
              )
          #print("radec polar")
          #print(f"{c_radec.ra}, {c_radec.dec}, {c_radec.distance}")
          c_radec_car = c_radec.cartesian
          #print("radec cartesian")
          #print(f"{c_radec_car.x}, {c_radec_car.y}, {c_radec_car.z}")

          c_g_eclip_mean = c_radec.transform_to(
              astropy.coordinates.builtin_frames.GeocentricMeanEcliptic)
          #print("g_eclip converted from radec polar")
          #print(f"{c_g_eclip_mean.lon}, {c_g_eclip_mean.lat}, {c_g_eclip_mean.distance}")
          c_g_eclip_mean = c_g_eclip_mean.cartesian
          x_earth = -c_g_eclip_mean.x.value
          y_earth = -c_g_eclip_mean.y.value
          z_earth = -c_g_eclip_mean.z.value

          x_sun_list.append(x_sun)
          y_sun_list.append(y_sun)
          z_sun_list.append(z_sun)
          x_earth_list.append(x_earth)
          y_earth_list.append(y_earth)
          z_earth_list.append(z_earth)

    dic_temp = dict(
        utc=utc_list, ra=ra_list, dec=dec_list,
        hecl_lon=hecl_lon, hecl_lat=hecl_lat, r=r_list, delta=delta_list,
        jd=jd_list, mjd=mjd_list,
        x_sun=x_sun_list, y_sun=y_sun_list, z_sun=z_sun_list,
        x_earth=x_earth_list, y_earth=y_earth_list, z_earth=z_earth_list
        )
    df_temp = pd.DataFrame(dic_temp.values(), index=dic_temp.keys()).T


    # function to predict various values when observation
    f_sun_x = interpolate.interp1d(
        df_temp["jd"], df_temp["x_sun"], kind='linear', fill_value='extrapolate')
    f_sun_y = interpolate.interp1d(
        df_temp["jd"], df_temp["y_sun"], kind='linear', fill_value='extrapolate')
    f_sun_z = interpolate.interp1d(
        df_temp["jd"], df_temp["z_sun"], kind='linear', fill_value='extrapolate')
    f_earth_x = interpolate.interp1d(
        df_temp["jd"], df_temp["x_earth"], kind='linear', fill_value='extrapolate')
    f_earth_y = interpolate.interp1d(
        df_temp["jd"], df_temp["y_earth"], kind='linear', fill_value='extrapolate')
    f_earth_z = interpolate.interp1d(
        df_temp["jd"], df_temp["z_earth"], kind='linear', fill_value='extrapolate')
    
    # Read photometric csv
    df = pd.read_csv(lc, sep=" ")
    col = df.columns.tolist()

    # predict x, y, z of the Sun and the Earth 
    df["x_sun"]   = df[key_jd].map(f_sun_x)
    df["y_sun"]   = df[key_jd].map(f_sun_y)
    df["z_sun"]   = df[key_jd].map(f_sun_z)
    df["x_earth"] = df[key_jd].map(f_earth_x)
    df["y_earth"] = df[key_jd].map(f_earth_y)
    df["z_earth"] = df[key_jd].map(f_earth_z)
    return df


def save4inv(result, absflux, random, key_jd, key_flux, key_fluxerr, out):
    """
    Save formatted text for convex inversion.

    Parameters
    ----------
    result : array of pandas.DataFrame
        array of formatted DataFrame
    absflux : bool
        whether flux is absolute
    random : bool
        whether shuffle data using fluxerr
    key_jd : str
        keyword of time in JD
    key_flux : str
        keyword of flux
    key_fluxerr : str
        keyword of flux uncertainty
    out : str
        output filename

    Return
    ------
    df_format : pandas.DataFrame

    """
    N = len(result)
    with open(out, "w") as f:
        # Number of lightcurves
        f.write(f"{N}\n")
        for n in range(N):
            # n-th lightcurve
            df = result[n]
            if random:
                # Monte Carlo 
                flux_mc = np.random.normal(
                    df[key_flux], df[key_fluxerr], len(df))
                df["flux1"] = flux_mc
            else:
                df["flux1"] = df[key_flux]

            # Absolute flux 1 (calibrated)
            if absflux:
                f.write(f"{len(df)} 1\n")
            # Relative flux 0
            else:
                f.write(f"{len(df)} 0\n")

            for idx,row in df.iterrows():
                if absflux:
                    flux = row["flux1"]
                else:
                    flux = row["flux1"]/np.median(df["flux1"])
                jd = row[key_jd]
                f.write(
                    f"{jd} {flux} "
                    f"{row['x_sun']} {row['y_sun']} {row['z_sun']} "
                    f"{row['x_earth']} {row['y_earth']} {row['z_earth']} "
                    f"\n")


def do_conv(lam, beta, lc):
    """
    Do convex inversion.

    Parameters
    ----------
    lam : float
        longitude
    beta : float
        latitude
    lc : str
        lightcurve
    """
    l = lam
    b = beta
    print(f"  (lam, beta) = ({l:.1f}, {b:.1f})")

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
