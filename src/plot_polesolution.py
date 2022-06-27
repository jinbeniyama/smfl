#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot convex inversion results.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  


if __name__ == "__main__":
    parser = ap(description="Plot result of convex inversion.")
    parser.add_argument(
        "--Nrow", type=int, default=1, 
        help="number of row includes chi square")
    parser.add_argument(
        "--Nlam", type=int, default=1, 
        help="number of ecliptic longitude")
    parser.add_argument(
        "--Nbeta", type=int, default=1, 
        help="number of ecliptic latitude")
    parser.add_argument(
        "--norm", action="store_true", default=False,
        help="Normalize chi2 by the minimum")
    parser.add_argument(
        "--outtype", type=str, default="png",
        help="format of output figures")
    args = parser.parse_args()
 
    # Color map
    cm = "inferno"

    # Number of data
    Nlam = args.Nlam
    Nbeta = args.Nbeta

    # Normalization
    norm = args.norm

    lam = np.linspace(0, 360, Nlam)
    beta = np.linspace(-90, 90, Nbeta)

    # Make grid
    xx, yy = np.meshgrid(lam, beta)
    data = np.c_[xx.ravel(), yy.ravel()]

    # List for chi2 values
    Z = []
    for x in data:
        l, b = x[0], x[1]
        res = f"res_ci_{int(l)}_{int(b)}"
        with open(res, "r") as f:
            lines = f.readlines()
            line = lines[args.Nrow-1:args.Nrow]
            # Extract chi2
            line = line[0].split("chi2")[-1]
            line = line.split("dev")[0]
            chi2 = line.strip(" ")

        Z.append(float(chi2))
    if norm:
        Z = [x/np.min(Z) for x in Z]
    Z = np.array(Z)
    Z = Z.reshape(xx.shape)

    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_axes([0.15, 0.15, 0.63, 0.8])
    ax1.set_xlabel(r"$\lambda$ [deg]")
    ax1.set_ylabel(r"$\beta$ [deg]")

    im = ax1.contourf(xx, yy, Z, zorder=-1, cmap=cm)
    cax = fig.add_axes([0.80, 0.15, 0.03, 0.8])
    cbar = fig.colorbar(
        im,  cax=cax, orientation='vertical', label=r"$\chi^2$")
    out = f"pole_solution_Nlam{Nlam}_Nbeta{Nbeta}_norm{norm}.{args.outtype}"
    plt.savefig(out)
