#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot sidereal period vs. chi square.

If len(out_period_scan) > 1 (i.e., N_mc > 0),
the uncertainty of rotation period is given by its standard deviation.
If len(out_period_scan) == 1 (i.e., N_mc = 1),
the uncertainty of rotation period is estimated using chi-square values.
The same approach in Fatka+2025, A&A, is used in the code.
This is 3-sigma-"like" uncertainty.
"""
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  

from calcerror import round_error
from smfl import plotstyle


def calc_confidence_chi2(p_list, chi2_list, dof):
    """Estimate rotation period and its error.

    Parameters
    ----------
    p_list : array-like
        list of rotation periods
    chi2_list : array-like
        corresponding chi2 values
    dof : int
        degrees of freedom
    
    Returns
    -------
    P_cand : array-like
        rotation period with chi2 smaller than chi2_3sigma
    chi2_cand : float
        chi-squared values corresponding to P_cand
    chi2_3sigma : float
        chi-squared boundary with 3-sigma confidence level
    """
    
    # Normalize chi2 with the best value
    chi2_min = np.min(chi2_list)
    chi2_norm_list = [x/chi2_min for x in chi2_list]

    # 3-sigma like boundary
    # confidence level of 0.9973
    # Since chi2 min is normalized to the unity above,
    # this boundary corresponds to those of Polishook 2014, Icarus, 241, 79 
    # and Cambioni+2021, Nature.       
    chi2_min_norm = np.min(chi2_norm_list)
    chi2_3sigma = chi2_min_norm + 3*(2/dof)**0.5

    # Search periods with chi2 values smaller than chi2_3sigma
    P_cand, chi2_cand = [], []
    for p, c in zip(p_list, chi2_norm_list):
        if c < chi2_3sigma:
            P_cand.append(p)
            chi2_cand.append(c)
    return P_cand, chi2_cand, chi2_3sigma


def plot_chi2_rotP(out_period_scan, obj, dof, Psec=False, out="chi2_rotP.jpg"):
    """Plot chi2 vs. rotP.

    Parameters
    ----------
    out_period_scan : array-like
        output of period_scan (N=1)
    obj : str
        object name
    dof : int
        degrees of freedom
    Psec : bool, optional
        whether the unit of rotation period is second
    out : str, bool 
        output file name
    """

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.80-0.05])
    # Margin for object name in title
    if Psec:
        ax.set_xlabel("Sidereal period [s]")
    else:
        ax.set_xlabel("Sidereal period [h]")
    ylabel = "Reduced $\chi^2$ (min is normalized to 1)"
    ax.set_ylabel(ylabel)
  
    offset = 0
    p_best_list = []
    p_list, rms_list, chi2_list, iter_list, dark_list = [], [], [], [], []
    with open (out_period_scan[0], "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.replace("\n", "")
            line = line.split(" ")
            # Remove ""
            line = [float(x) for x in line if x != ""]
            assert len(line)==5, "Check the input."
            if Psec:
                line[0] = line[0]*3600.
            p_list.append(line[0])
            rms_list.append(line[1])
            chi2_list.append(line[2])
            iter_list.append(line[3])
            dark_list.append(line[4])
    
    # Order by rotP
    zip_lists = zip(p_list, chi2_list)
    zip_sort = sorted(zip_lists)
    p_list, chi2_list = zip(*zip_sort)
    # Normalize chi2 with the best value
    chi2_min = np.min(chi2_list)
    chi2_norm_list = [x/chi2_min for x in chi2_list]

    col, _ = plotstyle(0)
    ax.scatter(p_list, chi2_norm_list, s=30, color=col, marker="x")
    
    # Estimate uncertainty of rotation period 
    # The same approach with Fatka+2025
    # Obtain normalized chi2 list and 3-sigma boundary
    # If there are more than two candidates,
    # the analysis would be complex.
    P_cand, chi2_cand, chi2_3sigma = calc_confidence_chi2(p_list, chi2_list, dof)
    # Plot candidates
    for p, c in zip(P_cand, chi2_cand):
        ax.scatter(p, c,  color="blue", s=100, fc="None", marker="o")

    ax.set_title(f"{obj} (dof={dof})")
    plt.savefig(out)
    plt.close()


def plot_chi2_rotP_MC(out_period_scan, obj, Psec=False, out="chi2_rotP_MC.jpg"):
    """Plot chi2 vs. rotP with Monte-Carlo approach.

    Parameters
    ----------
    out_period_scan : array-like
        output of period_scan
    obj : str
        object name
    Psec : bool, optional
        whether the unit of rotation period is second
    out : str, bool 
        output file name
    """

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.80-0.05])
    # Margin for object name in title
    if Psec:
        ax.set_xlabel("Sidereal period [s]")
    else:
        ax.set_xlabel("Sidereal period [h]")
    ylabel = "Reduced $\chi^2$ (+offset)"
    ax.set_ylabel(ylabel)
  
    offset = 0
    p_best_list = []
    for n in range(N_mc):
        p_list, rms_list, chi2_list, iter_list, dark_list = [], [], [], [], []
        with open (out_period_scan[n], "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.replace("\n", "")
                line = line.split(" ")
                # Remove ""
                line = [float(x) for x in line if x != ""]
                assert len(line)==5, "Check the input."
                if Psec:
                    line[0] = line[0]*3600.
                p_list.append(line[0])
                rms_list.append(line[1])
                chi2_list.append(line[2])
                iter_list.append(line[3])
                dark_list.append(line[4])
        
        # Order by rotP
        zip_lists = zip(p_list, chi2_list)
        zip_sort = sorted(zip_lists)
        p_list, chi2_list = zip(*zip_sort)

        col, _ = plotstyle(n)
        # Shift to match the first result
        if n!=0:
            offset = min0 - np.min(chi2_list)
        chi2_list = [x+offset for x in chi2_list]
        if N_mc == 1:
            ax.scatter(p_list, chi2_list, s=30, color=col, marker="x")
        else:
            ax.plot(p_list, chi2_list, lw=1, color=col)
        
        # Save the minimum value of the first try
        if n == 0:
            min0 = np.min(chi2_list)
  
        # Save period
        idx_min = chi2_list.index(min(chi2_list))
        p_best_list.append(p_list[idx_min])

    # Estimate uncertainty of rotation period 
    # This is qualitativly different from the conventinal one..
    rotP, rotPerr = np.median(p_best_list), np.std(p_best_list)

    # Plot histograms
    print(f"p_mean, p_std = {p_mean}, {p_std}")
    rotP_str,rotPerr_str = round_error(rotP, rotPerr)
    print(f"p_mean, p_std = {rotP_str}, {rotPerr_str}")
    if args.sec:
        Ppart = f"P=${p_mean_str}\pm{p_std_str} s$"
    else:
        Ppart = f"P=${p_mean_str}\pm{p_std_str} h$"
    
    # Plot N_mc and estimated rotP
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x = xmin + (xmax-xmin)*0.05
    y = ymin + (ymax-ymin)*0.1
    ax.text(x, y, f"N={N_mc}\n{Ppart}", fontsize=25)
  
    ymin, ymax = ax.get_ylim()
    ax.vlines(p_mean, ymin, ymax, color="black", ls="solid")
    ax.vlines(p_mean-p_std, ymin, ymax, color="black", ls="dashed")
    ax.vlines(p_mean+p_std, ymin, ymax, color="black", ls="dashed")
    ax.set_ylim([ymin, ymax])
    
    ax.set_title(f"{obj}")
    plt.savefig(out)
    plt.close()




if __name__ == "__main__":
    parser = ap(description="Plot period vs. chi2.")
    parser.add_argument(
        "obj", type=str, 
        help="object name")
    parser.add_argument(
        "dof", type=int, 
        help="Degree of freedom")
    parser.add_argument(
        "out_period_scan", nargs="*", 
        help="output file of period_scan")
    parser.add_argument(
        "--sec", action="store_true", default=False,
        help="Handle period in sec")
    parser.add_argument(
        "--out", type=str, default="ps_result.png",
        help="Output file name")
    args = parser.parse_args()
    
    # Number of trials with Monte Carlo technique
    N_mc = len(args.out_period_scan)

    if N_mc > 1:
        plot_chi2_rotP_MC(args.out_period_scan, args.obj, args.sec, out=args.out)
    else:
        plot_chi2_rotP(args.out_period_scan, args.obj, args.dof, args.sec, out=args.out)
