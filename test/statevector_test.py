#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test of how to get state vectors (Earth, Sun).

Method 1 and 2 give diffrent results.
Method 1 seems better.
"""
from argparse import ArgumentParser as ap
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy


def state_vec(target, jd, code, method):
    
    if method == 1:
        S = Horizons(location="500@10", id=target, epochs=jd)
        vec = S.vectors(refplane="ecliptic")
        # Vector from the Sun to asteroid (Sun -> ast)
        x_S, y_S, z_S = vec["x"][0], vec["y"][0], vec["z"][0]
        # Vector from asteroid to the Sun (ast -> Sun)
        x_S, y_S, z_S = -x_S, -y_S, -z_S

        # Location of code means observer center
        E = Horizons(location=code, id=target, epochs=jd)
        vec = E.vectors(refplane="ecliptic")
        eph = E.ephemerides()

        # Vector from the Earth to asteroid (Earth -> ast)
        x_E, y_E, z_E = vec["x"][0], vec["y"][0], vec["z"][0]
        # Vector from asteroid to the Earth (ast -> Earth)
        # Vector from asteroid to the Earth (ast -> Earth)
        x_E, y_E, z_E = -x_E, -y_E, -z_E


    # Old method by JB
    if method == 2:
        ast = Horizons(location=code, id=target, epochs=jd)
        eph = ast.ephemerides()
        lam, beta, r   = eph["EclLon"], eph["EclLat"], eph["r"]
        ra, dec, delta = eph["RA"], eph["DEC"], eph["delta"]

        ## use cartisian method
        c_eclip = SkyCoord(
            lam, beta, r, frame="heliocentricmeanecliptic", 
            unit=(u.deg, u.deg, u.au)
            )
        c_eclip = c_eclip.cartesian

        x_S = -c_eclip.x.value[0]
        y_S = -c_eclip.y.value[0]
        z_S = -c_eclip.z.value[0]

        # not icrs (), but gcrs(Geocentric Celestial Reference System)
        c_radec = SkyCoord(
            ra=ra, dec=dec, distance=delta, frame="gcrs", unit=(u.hourangle, u.deg, u.au),
        )
        c_radec_car = c_radec.cartesian
        #print("radec cartesian")
        #print(f"{c_radec_car.x}, {c_radec_car.y}, {c_radec_car.z}")

        c_g_eclip_mean = c_radec.transform_to(
            astropy.coordinates.builtin_frames.GeocentricMeanEcliptic)
        c_g_eclip_mean = c_g_eclip_mean.cartesian
        x_E = -c_g_eclip_mean.x.value[0]
        y_E = -c_g_eclip_mean.y.value[0]
        z_E = -c_g_eclip_mean.z.value[0]

    return x_E, y_E, z_E, x_S, y_S, z_S



if __name__ == "__main__":
    parser = ap(description="Test of how to get state vectors.")
    args = parser.parse_args()
    
    # Time before light-time corrected
    t_jd0 = 2458521.53767
    # Time after light-time corrected
    t_jd1 = 2458521.53670
    # Light-traveling time is ~ 83 s
    target = "2015 BY310"

    code_list = ["500", "W74"]
    code_list = ["500"]

    # Let's get state vectors
    print("")
    print("Without light-time correction")
    for co in code_list:
        for n_method in range(1, 3, 1):
            x_E, y_E, z_E, x_S, y_S, z_S = state_vec(target, t_jd0, co, method=n_method)
            print(f"  Method {n_method} at {co} (code)")
            print(f"    x_S, y_S, z_S = {x_S}, {y_S}, {z_S}")
            print(f"    x_E, y_E, z_E = {x_E}, {y_E}, {z_E}")

    print("")
    print("")
    
    print("With light-time correction")
    for co in code_list:
        for n_method in range(1, 3, 1):
            x_E, y_E, z_E, x_S, y_S, z_S = state_vec(target, t_jd1, co, method=n_method)
            print(f"  Method {n_method} at {co} (code)")
            print(f"    x_S, y_S, z_S = {x_S}, {y_S}, {z_S}")
            print(f"    x_E, y_E, z_E = {x_E}, {y_E}, {z_E}")
    


