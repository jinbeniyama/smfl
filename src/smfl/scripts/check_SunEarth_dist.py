#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot light curve with modeled one.

To be updated : Calculate and plot phase from jd.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  


if __name__ == "__main__":
    parser = ap(description="Check data handling.")
    parser.add_argument(
        "csv1", help="jd flux x y z (Sun) x y z (Earth)")
    parser.add_argument(
        "csv2", help="jd flux x y z (Sun) x y z (Earth)")
    parser.add_argument(
        "Nline", type=int,
        help="Number of lines to be used.")
    args = parser.parse_args()
    
    
    # Read data 1
    jd1_list = []
    xs1_list = []
    ys1_list = []
    zs1_list = []
    xe1_list = []
    ye1_list = []
    ze1_list = []
    with open(args.csv1, "r") as f:
        lines = f.readlines()
        lines = lines[2:2+args.Nline]
        for l in lines:
            l = l.split(" ")
            l = [float(x) for x in l if x!="\n"]
            jd1_list.append(l[0])
            xs1_list.append(l[2])
            ys1_list.append(l[3])
            zs1_list.append(l[4])
            xe1_list.append(l[5])
            ye1_list.append(l[6])
            ze1_list.append(l[7])

    # Read data 2
    jd2_list = []
    xs2_list = []
    ys2_list = []
    zs2_list = []
    xe2_list = []
    ye2_list = []
    ze2_list = []
    with open(args.csv2, "r") as f:
        lines = f.readlines()
        lines = lines[2:2+args.Nline]
        for l in lines:
            l = l.split(" ")
            l = [float(x) for x in l if (x!="\n") and (x!="")]
            jd2_list.append(l[0])
            xs2_list.append(l[2])
            ys2_list.append(l[3])
            zs2_list.append(l[4])
            xe2_list.append(l[5])
            ye2_list.append(l[6])
            ze2_list.append(l[7])
    print(xs1_list)
    print(xs2_list)
    fig = plt.figure(figsize=(10, 16))
    ax1_1   = fig.add_axes([0.15, 0.8, 0.3, 0.15])
    ax1_1_d = fig.add_axes([0.15, 0.73, 0.3, 0.05])
    ax1_2   = fig.add_axes([0.55, 0.8, 0.3, 0.15])
    ax1_2_d = fig.add_axes([0.55, 0.73, 0.3, 0.05])

    ax2_1   = fig.add_axes([0.15, 0.5, 0.3, 0.15])
    ax2_1_d = fig.add_axes([0.15, 0.43, 0.3, 0.05])
    ax2_2 = fig.add_axes([0.55, 0.5, 0.3, 0.15])
    ax2_2_d = fig.add_axes([0.55, 0.43, 0.3, 0.05])

    ax3_1   = fig.add_axes([0.15, 0.2, 0.3, 0.15])
    ax3_1_d = fig.add_axes([0.15, 0.13, 0.3, 0.05])
    ax3_2   = fig.add_axes([0.55, 0.2, 0.3, 0.15])
    ax3_2_d = fig.add_axes([0.55, 0.13, 0.3, 0.05])
    

    ax3_1.set_xlabel("JD [day]")
    ax3_2.set_xlabel("JD [day]")
    ax1_1.set_ylabel("Distance [au]")
    ax2_1.set_ylabel("Distance [au]")
    ax3_1.set_ylabel("Distance [au]")

    ax1_1.scatter(
        jd1_list, xs1_list, color="red", marker="o", 
        label="Sun x 1", ls="solid")
    ax1_1.scatter(
        jd2_list, xs2_list, color="red", marker="x", 
        label="Sun x 2", ls="dashed")
    ax2_1.scatter(
        jd1_list, ys1_list, color="blue", marker="o", 
        label="Sun y 1", ls="solid")
    ax2_1.scatter(
        jd2_list, ys2_list, color="blue", marker="x", 
        label="Sun y 2", ls="dashed")
    ax3_1.scatter(
        jd1_list, zs1_list, color="green", marker="o", 
        label="Sun z 1", ls="solid")
    ax3_1.scatter(
        jd2_list, zs2_list, color="green", marker="x", 
        label="Sun z 2", ls="dashed")

    ax1_2.scatter(
        jd1_list, xe1_list, color="red", marker="o", 
        label="Earth x 1", ls="solid")
    ax1_2.scatter(
        jd2_list, xe2_list, color="red", marker="x", 
        label="Earth x 2", ls="dashed")
    ax2_2.scatter(
        jd1_list, ye1_list, color="blue", marker="o", 
        label="Earth y 1", ls="solid")
    ax2_2.scatter(
        jd2_list, ye2_list, color="blue", marker="x", 
        label="Earth y 2", ls="dashed")
    ax3_2.scatter(
        jd1_list, ze1_list, color="green", marker="o", 
        label="Earth z 1", ls="solid")
    ax3_2.scatter(
        jd2_list, ze2_list, color="green", marker="x", 
        label="Earth z 2", ls="dashed")


    # Plot dirrerence
    xs_diff = [(x-y)*1e8 for x,y in zip(xs1_list, xs2_list)]
    ys_diff = [(x-y)*1e8 for x,y in zip(ys1_list, ys2_list)]
    zs_diff = [(x-y)*1e8 for x,y in zip(zs1_list, zs2_list)]
    xe_diff = [(x-y)*1e8 for x,y in zip(xe1_list, xe2_list)]
    ye_diff = [(x-y)*1e8 for x,y in zip(ye1_list, ye2_list)]
    ze_diff = [(x-y)*1e8 for x,y in zip(ze1_list, ze2_list)]

    ax1_1_d.scatter(
        jd1_list, xs_diff,marker="o", 
        label="Sun x diff", ls="solid")
    ax1_2_d.scatter(
        jd1_list, xe_diff,marker="o", 
        label="Earth x diff", ls="solid")
    ax2_1_d.scatter(
        jd1_list, ys_diff,marker="o", 
        label="Sun y diff", ls="solid")
    ax2_2_d.scatter(
        jd1_list, ye_diff,marker="o", 
        label="Earth y diff", ls="solid")
    ax3_1_d.scatter(
        jd1_list, zs_diff,marker="o", 
        label="Sun z diff", ls="solid")
    ax3_2_d.scatter(
        jd1_list, ze_diff,marker="o", 
        label="Earth z diff", ls="solid")

    for ax in fig.axes:
        ax.legend()

    out = f"check_SunEarth.png"
    plt.savefig(out)
    plt.close()
