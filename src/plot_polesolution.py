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
        "--N_chi2", type=int, default=1, 
        help="number of row includes chi square")
    parser.add_argument(
        "--N_P", type=int, default=1, 
        help="number of row includes rotP")
    parser.add_argument(
        "--N_dfac", type=int, default=1, 
        help="number of row includes dark facet area")
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
    
    N_chi2 = args.N_chi2
    N_P = args.N_P
    N_dfac = args.N_dfac

    # Normalization
    norm = args.norm

    lam = np.linspace(0, 360, Nlam)
    beta = np.linspace(-90, 90, Nbeta)

    # Make grid
    xx, yy = np.meshgrid(lam, beta)
    data = np.c_[xx.ravel(), yy.ravel()]

    # List for chi2 values
    chi2_list = []
    rotP_list = []
    dfac_list = []
    for x in data:
        l, b = x[0], x[1]
        res = f"res_ci_{int(l)}_{int(b)}"
        with open(res, "r") as f:
            lines = f.readlines()

            # Extract chi2
            line = lines[N_chi2-1:N_chi2]
            line = line[0].split("chi2")[-1]
            line = line.split("dev")[0]
            chi2 = line.strip(" ")

            # Extract rotP
            line = lines[N_P-1:N_P]
            line = line[0].split(" ")
            rotP = line[-1].strip("\n")
            # hour to sec
            rotP = float(rotP)*3600.

            # Extract dfac (dark facet area)
            line = lines[N_dfac-1:N_dfac]
            line = line[0].split(" ")
            dfac = line[-1].strip("%\n")

            chi2_list.append(float(chi2))
            rotP_list.append(float(rotP))
            dfac_list.append(float(dfac))
    if norm:
        Z_chi2 = [x/np.min(chi2_list) for x in chi2_list]
        # percentage from the minima
        Z_chi2 = [(x-1)*100 for x in Z_chi2]

        Z_chi2 = np.array(Z_chi2)
        #Z_rotP = np.array(Z_rotP)
        label_chi2 = r"$\chi^2$ [%]"
    else:
        Z_chi2 = np.array(chi2_list)
        label_chi2 = r"$\chi^2$"

    Z_rotP = np.array(rotP_list)
    Z_dfac = np.array(dfac_list)

    Z_chi2 = Z_chi2.reshape(xx.shape)
    Z_rotP = Z_rotP.reshape(xx.shape)
    Z_dfac = Z_dfac.reshape(xx.shape)

    fig = plt.figure(figsize=(8, 14))

    ax1 = fig.add_axes([0.15, 0.72, 0.63, 0.25])
    cax1 = fig.add_axes([0.80, 0.72, 0.03, 0.25])

    ax2 = fig.add_axes([0.15, 0.40, 0.63, 0.25])
    cax2 = fig.add_axes([0.80, 0.40, 0.03, 0.25])

    ax3 = fig.add_axes([0.15, 0.08, 0.63, 0.25])
    cax3 = fig.add_axes([0.80, 0.08, 0.03, 0.25])

    # ToDo : Avoid "+1e5" etc. 
    import matplotlib
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    cax2.yaxis.set_major_formatter(y_formatter)
    cax2.xaxis.set_major_formatter(y_formatter)
    
    ax1.set_xlabel(r"$\lambda$ [deg]")
    ax2.set_xlabel(r"$\lambda$ [deg]")
    ax3.set_xlabel(r"$\lambda$ [deg]")
    ax1.set_ylabel(r"$\beta$ [deg]")
    ax2.set_ylabel(r"$\beta$ [deg]")
    ax3.set_ylabel(r"$\beta$ [deg]")

    im = ax1.contourf(xx, yy, Z_chi2, zorder=-1, cmap=cm)
    cbar = fig.colorbar(
        im,  cax=cax1, orientation='vertical', label=label_chi2)

    im = ax2.contourf(xx, yy, Z_rotP, zorder=-1, cmap=cm)
    cbar = fig.colorbar(
        im,  cax=cax2, orientation='vertical', label=r"$Period [s]$")

    im = ax3.contourf(xx, yy, Z_dfac, zorder=-1, cmap=cm)
    cbar = fig.colorbar(
        im,  cax=cax3, orientation='vertical', label="Dark facet area [%]")

    out = f"pole_solution_Nlam{Nlam}_Nbeta{Nbeta}_norm{norm}.{args.outtype}"
    plt.savefig(out)
