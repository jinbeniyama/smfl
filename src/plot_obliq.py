#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot obliquity.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
from astroquery.jplhorizons import Horizons


def calc_obliq(lam, beta, incl, LoA):
    """
    Calculate obliquity from ecliptic coordinates.

    Parameter
    ---------
    lam, beta : float
        ecliptic coordinates of poles in deg
    incl : float
        orbital inclination in deg
    LoA : float
        longitude of ascending node in deg

    Return
    ------
    obliq : float
        obliquity in deg
    """

    # Orbital elements
    sin_i = np.sin(np.deg2rad(incl))
    cos_i = np.cos(np.deg2rad(incl))
    sin_Omega = np.sin(np.deg2rad(LoA))
    cos_Omega = np.cos(np.deg2rad(LoA))

    # Orbit pole vector
    x_orb = sin_i*sin_Omega
    y_orb = -cos_Omega*sin_i
    z_orb = cos_i

    # Ecliptic coordinates
    sin_lambda = np.sin(np.deg2rad(lam))
    cos_lambda = np.cos(np.deg2rad(lam))
    sin_beta = np.sin(np.deg2rad(beta))
    cos_beta = np.cos(np.deg2rad(beta))

    # Asteroid spin pole vector
    x_spin = cos_beta*cos_lambda
    y_spin = cos_beta*sin_lambda
    z_spin = sin_beta

    # Inner product equals cos(obliq)
    cos_obliq = x_spin*x_orb + y_spin*y_orb + z_spin*z_orb
    # In degree
    obliq = np.rad2deg(np.arccos(cos_obliq))
    return obliq


if __name__ == "__main__":
    parser = ap(description="Plot obliquity map.")
    parser.add_argument(
        "obj", type=str,  
        help="Object name")
    parser.add_argument(
        "--lam", type=float, default=10, 
        help="pole longitude in deg")
    parser.add_argument(
        "--beta", type=float, default=10, 
        help="pole latitude in deg")
    parser.add_argument(
        "--outtype", type=str, default="jpg",
        help="format of output figures")
    args = parser.parse_args()
 
    # Color map
    cm = "RdYlBu"

    # Number of data
    Nlam = 360
    Nbeta = 180
    
    lam = np.linspace(0, 360, Nlam)
    beta = np.linspace(-90, 90, Nbeta)
    # Make grid
    xx, yy = np.meshgrid(lam, beta)
    data = np.c_[xx.ravel(), yy.ravel()]

    # Obtain inclination and LoA to calculate obliquity
    ast = Horizons(id=args.obj, location="500@10", epochs=2450951.5)
    el = ast.elements()
    incl = el["incl"][0]
    LoA = el["Omega"][0]

    obliq_list = []
    for x in data:
        l, b = x[0], x[1]
        # Calculate obliquity
        obliq = calc_obliq(l, b, incl, LoA)
        obliq_list.append(obliq)

    lam = np.linspace(-180, 180, Nlam)
    beta = np.linspace(-90, 90, Nbeta)

    # Make grid
    xx, yy = np.meshgrid(lam, beta)
    data = np.c_[xx.ravel(), yy.ravel()]
    Z_obliq = np.array(obliq_list)
    Z_obliq = Z_obliq.reshape(xx.shape)

    fig = plt.figure(figsize=(8, 6))
    
    ax1 = fig.add_axes([0.15, 0.20, 0.63, 0.70])
    cax1 = fig.add_axes([0.80, 0.20, 0.03, 0.70])

    ax1.set_xlabel(r"$\lambda$ [deg]")
    ax1.set_ylabel(r"$\beta$ [deg]")
    
    label_obliq = "Obliquity [deg]"
    levels = np.arange(0, 181, 10)
    im = ax1.contourf(xx, yy, Z_obliq, levels, cmap=cm)
    cbar = fig.colorbar(
        im,  cax=cax1, orientation='vertical', label=label_obliq)

    # Plot pole
    if (args.lam) and (args.beta):
        label_pole = r"Pole ($\lambda, \beta$) = " + f"({args.lam}, {args.beta})" 
        ax1.scatter(
            args.lam, args.beta, marker="*", color="green", s=400, ec="black", 
            label=label_pole)
    ax1.legend().get_frame().set_alpha(1.0)
    out = f"obliq.jpg"
    plt.savefig(out)
