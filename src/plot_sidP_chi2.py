#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot sidereal period vs. chi square.
"""
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
import datetime

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
      "--outtype", type=str, default="png",
      help="format of output figures")
    args = parser.parse_args()
    
    # Number of trials with Monte Carlo technique
    N_mc = len(args.out_period_scan)
  
  
    fig = plt.figure()
    ax = fig.add_axes([0.2, 0.15, 0.75, 0.65])
    ax_p_hist = fig.add_axes([0.2, 0.80, 0.750, 0.1])
    ax.set_xlabel("Sidereal period [h]")
    if N_mc == 0:
        ylabel = "$\chi^2$"
    else:
        ylabel = "$\chi^2$ (shifted)"
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
                p_list.append(line[0])
                rms_list.append(line[1])
                chi2_list.append(line[2])
                iter_list.append(line[3])
                dark_list.append(line[4])
  
        col, mark = plotstyle(n)
        print(f"  col, mark = {col}, {mark}")
        # Shift to match the first result
        if n!=0:
            offset = min0 - np.min(chi2_list)
        chi2_list = [x+offset for x in chi2_list]
        ax.scatter(p_list, chi2_list, s=7, color=col, marker=mark)
        
        # Save the minimum value of the first try
        if n==0:
            min0 = np.min(chi2_list)
  
        # Save period
        idx_min = chi2_list.index(min(chi2_list))
        p_best_list.append(p_list[idx_min])
  
  
    # Plot histograms
    p_mean, p_std = np.mean(p_best_list), np.std(p_best_list)
    print(f"p_mean, p_std = {p_mean}, {p_std}")
    p_mean_str, p_std_str = round_error(p_mean, p_std)
    print(f"p_mean, p_std = {p_mean_str}, {p_std_str}")
    label = f"P=${p_mean_str}\pm{p_std_str} h$"
    ax_p_hist.hist(
      p_best_list, color="black", label=label)
    ax_p_ymax = ax_p_hist.get_ylim()[1]
    ax_p_hist.vlines(p_mean, 0, ax_p_ymax, color="black", ls="solid")
    ax_p_hist.vlines(p_mean-p_std, 0, ax_p_ymax, color="black", ls="dashed")
    ax_p_hist.vlines(p_mean+p_std, 0, ax_p_ymax, color="black", ls="dashed")
    ax_p_hist.axes.xaxis.set_visible(False)
    pmin, pmax = ax.get_xlim()
    ax_p_hist.set_xlim([pmin, pmax])
  
    ax.set_title(f"{args.obj}")
    ax.legend()
    ax_p_hist.legend()
    out = f"{args.obj}_sidP_chi2_N{N_mc}.{args.outtype}"
    plt.savefig(out)
    plt.close()
