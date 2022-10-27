#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot convex inversion results.
Uncertainty is estimated in the same manner as the previous papers
such as Vokrouhlicky et al. (2011, 2017) and Durech et al. (2018).
The lightcurve is needed as an argument to calculate the number of observations N.


TODO
----
To calculate the number of degree of freedom, nu = N - M,
this code assumes M = 100 (i.e., nu = N - 100).
This does not change the result for normal shape modeling (N >> 100),
but I personally feel the exact number (M and N) should be writen 
in the manuscript if you are writing a paper.


References
----------
Vokrouhlicky et al. 2011, AJ, 142, 159.
Vokrouhlicky et al. 2017, AJ, 153, 270.
Durech et al. 2018, AA, 609, A86.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib

from smfl import nobs_lc


def chi2_CI(N=2, M=1, nu=None, sigma=3):
    """
    Calculate confidense interval (CI) of chi2 distribution.
    
    Parameters
    ----------
    N : int
      number of observations
    M : int
      number of parameters in convex inv
    nu : int, optional 
      degree of freedom of chi2 distribution (N-M)
    sigma : int
      confidence level

    Return
    ------
    CI : float
      n-sigma confidence interval in percentage
    """
    if not nu:
        nu = N - M
    # For reduced chi2
    CI = sigma*np.sqrt(2/nu)
    ## For normal chi2
    #CI = sigma*np.sqrt(2*nu)
    return CI


if __name__ == "__main__":
    parser = ap(description="Plot result of convex inversion.")
    parser.add_argument(
        "--Nlam", type=int, default=1, 
        help="number of ecliptic longitude")
    parser.add_argument(
        "--Nbeta", type=int, default=1, 
        help="number of ecliptic latitude")
    parser.add_argument(
        "--nu", type=int, default=1000, 
        help="degree of freedom ")
    parser.add_argument(
        "--lc", type=str, default=None, 
        help="Formatted lightcurve to calculate the number of obs.")
    parser.add_argument(
        "--sigma", type=float, default=3, 
        help="CI level")
    parser.add_argument(
        "--norm", action="store_true", default=False,
        help="Normalize chi2 by the minimum")
    parser.add_argument(
        "--outtype", type=str, default="png",
        help="format of output figures")
    parser.add_argument(
        "--resdir", type=str, default="convex_result",
        help="Directory for results of convexinv")
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
    chi2_list = []
    rotP_list = []
    dfac_list = []
    for x in data:
        l, b = x[0], x[1]
        res = f"res_ci_{int(l)}_{int(b)}"
        res = os.path.join(args.resdir, res)
        with open(res, "r") as f:

            # Extract last 2 + 5 columns
            lines = f.readlines()[-7:]
            # === example ===
            # 50  chi2 4.710719  dev 0.055289  alambda 0.001000
            #
            # lambda, beta and period (hrs): 360.000000 -90.000000 0.014347
            # phase function parameters: 0.000000 0.100000 0.000000 
            # Lambert coefficient: 0.1
            # plus a dark facet with area 1.28%
            # === example ===

            # Extract chi2
            line = lines[0]
            line = line.split("chi2")[-1]
            line = line.split("dev")[0]
            chi2 = line.strip(" ")

            # Extract rotP
            line = lines[2]
            print(line)
            line = line.split(" ")
            rotP = line[-1].strip("\n")
            # hour to sec
            rotP = float(rotP)*3600.

            # Extract dfac (dark facet area)
            line = lines[5]
            line = line.split(" ")
            dfac = line[-1].strip("%\n")

            chi2_list.append(float(chi2))
            rotP_list.append(float(rotP))
            dfac_list.append(float(dfac))


    chi2_min = np.min(chi2_list)
    if norm:
        Z_chi2 = [x/chi2_min for x in chi2_list]
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
    
    # How to calculate nu ?
    # Assume nu = N - M ~ N (i.e., N >> M)
    if args.lc:
        N = nobs_lc(args.lc)
        nu = N
    else:
        nu = args.nu
    CI = chi2_CI(nu=nu, sigma=args.sigma)
    print(f"{args.sigma}-sigma CI is {CI:.2f}")
    # Where from?
    print(f"{args.sigma}-sigma range is chi2_min to chi2_min(1 + CI)")
    print(f"  -> chi2: {chi2_min:.2f} to {chi2_min*(1 + CI):.2f}")
    # Add CI line
    levels = [chi2_min*(1 + CI)]
    ax1.contour(xx, yy, Z_chi2, levels, linestyles=["dashed"], colors="white")

    im = ax2.contourf(xx, yy, Z_rotP, zorder=-1, cmap=cm)
    cbar = fig.colorbar(
        im,  cax=cax2, orientation='vertical', label=r"$Period [s]$")

    im = ax3.contourf(xx, yy, Z_dfac, zorder=-1, cmap=cm)
    cbar = fig.colorbar(
        im,  cax=cax3, orientation='vertical', label="Dark facet area [%]")

    out = f"pole_solution_Nlam{Nlam}_Nbeta{Nbeta}_norm{norm}.{args.outtype}"
    plt.savefig(out)
