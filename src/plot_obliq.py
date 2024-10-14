#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot obliquity.

TODO:
    Query incl. and LoA with object name.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib

from smfl import calc_obliq




if __name__ == "__main__":
    parser = ap(description="Plot obliquity map.")
    parser.add_argument(
        "incl", type=float, default=10, 
        help="Orbital inclination in deg")
    parser.add_argument(
        "LoA", type=float, default=10, 
        help="Longitude of ascending node in deg")
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

    obliq_list = []
    for x in data:
        l, b = x[0], x[1]
        # Calculate obliquity
        obliq = calc_obliq(l, b, args.incl, args.LoA)
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
