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


if __name__ == "__main__":
    parser = ap(description="Plot period vs. chi2.")
    parser.add_argument(
        "obj", type=str, help="object name")
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

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.80-0.05])
    # Margin for object name in title
    if args.sec:
        ax.set_xlabel("Sidereal period [s]")
    else:
        ax.set_xlabel("Sidereal period [h]")
    if N_mc == 0:
        ylabel = "$Reduced \chi^2$"
    else:
        ylabel = "Reduced $\chi^2$ (+offset)"
    ax.set_ylabel(ylabel)
  
    offset = 0
    p_best_list = []
    for n in range(N_mc):
        p_list, rms_list, chi2_list, iter_list, dark_list = [], [], [], [], []
        with open (args.out_period_scan[n], "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.replace("\n", "")
                line = line.split(" ")
                # Remove ""
                line = [float(x) for x in line if x != ""]
                assert len(line)==5, "Check the input."
                if args.sec:
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

  
  
    # Plot histograms
    p_mean, p_std = np.mean(p_best_list), np.std(p_best_list)
    print(f"p_mean, p_std = {p_mean}, {p_std}")
    p_mean_str, p_std_str = round_error(p_mean, p_std)
    print(f"p_mean, p_std = {p_mean_str}, {p_std_str}")
    if args.sec:
        Ppart = f"P=${p_mean_str}\pm{p_std_str} s$"
    else:
        Ppart = f"P=${p_mean_str}\pm{p_std_str} h$"

    # If N_mc > 1, plot the number of trials and period with uncertainty
    if N_mc > 1:
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
    
    ax.set_title(f"{args.obj}")
    plt.savefig(args.out)
    plt.close()
