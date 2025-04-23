#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot light curve with model curve.
"""
from argparse import ArgumentParser as ap
import pandas as pd
import matplotlib.pyplot as plt  

from smfl import figure4lc, mycolor


if __name__ == "__main__":
    parser = ap(description="Plot lightcurves")
    parser.add_argument(
        "obj", help="object name")
    parser.add_argument(
        "lc_phot", help="photometric light curve")
    parser.add_argument(
        "--lc_model", default=None,
        help="modeled light curve")
    parser.add_argument(
        "--deltat", type=str, default=1,
        help="Threshold to judge different observations in day")
    parser.add_argument(
        "--obj", type=str, default="Asteroid",
        help="object name")
    parser.add_argument(
        "--rotP", type=float, default=None,
        help="rotational period in sec by default")
    parser.add_argument(
        "--hour", action="store_true",
        help="rotational period in hour")
    parser.add_argument(
        "--JD0", default=2459702.5, type=float,
        help="time zero point in Juliand day (default is 2022-05-03)")
    parser.add_argument(
        "--ylim", type=float, nargs=2, default=False,
        help="Y min and max")
    args = parser.parse_args()

    
    # Time zero point in Julian day
    JD0 = args.JD0



    # Read photometric result =================================================
    jd_phot, flux_phot, n_lc = [], [], []
    # Index of the number of lightcurve
    idx_lc = 1
    with open(args.lc_phot, "r") as f:
        # Skip headers
        lines = f.readlines()[2:]
        for l in lines:
            l = l.split(" ")
            N_line = len(l)
            if N_line < 8:
                idx_lc += 1
                continue
            l = [float(x) for x in l if ((x!="\n") & (x!=""))]
            jd_phot.append(l[0])
            flux_phot.append(l[1])
            n_lc.append(idx_lc)
    df = pd.DataFrame(dict(jd=jd_phot, flux=flux_phot, n_lc=n_lc))
    # Number of lightcurves
    N_lc = idx_lc
    print(f"Number of lightcurves = {N_lc}")
    # Read photometric result =================================================


    # Read model curves (optional ) ===========================================
    # N_model should be N_phot
    if args.lc_model:
        with open(args.lc_model, "r") as f:
            lines = f.readlines()
            lc_model = [float(x.replace("\n", "")) for x in lines]
            assert len(df) == len(lc_model), "Check the input files."
            df["flux_model"] = lc_model
        str_model = "_wmodel"
    else:
        str_model = ""
    # Read model curves (optional ) ===========================================

   
    # TODO: update for arbtary lightucrves.
    #  Now optimized for 7(8) lightcurves !!
    if N_lc > 8:
        print("Only show 8 lightcurves.")


    # Plot lightcurves ========================================================
    # Create 4x2 figures lightcurves
    fig = plt.figure(figsize=(15, 21)) 
    ax1_1 = fig.add_axes([0.1, 0.80, 0.35, 0.18])
    ax1_2 = fig.add_axes([0.6, 0.80, 0.35, 0.18])
    ax2_1 = fig.add_axes([0.1, 0.55, 0.35, 0.18])
    ax2_2 = fig.add_axes([0.6, 0.55, 0.35, 0.18])
    ax3_1 = fig.add_axes([0.1, 0.30, 0.35, 0.18])
    ax3_2 = fig.add_axes([0.6, 0.30, 0.35, 0.18])
    ax4_1 = fig.add_axes([0.1, 0.05, 0.35, 0.18])
    ax4_2 = fig.add_axes([0.6, 0.05, 0.35, 0.18])

     
    for idx, ax in enumerate(
        [ax1_1, ax1_2, ax2_1, ax2_2, ax3_1, ax3_2, ax4_1, ax4_2]):
        ax.set_xlabel(f"JD-{JD0} [day]")
        ax.set_ylabel("Relative flux")
        df_temp = df[df["n_lc"]==idx+1]
        ax.scatter(
            df_temp["jd"]-JD0, df_temp["flux"], 
            color=mycolor[1], label=f"LC {idx+1} obs")

        if args.lc_model:
            ax.scatter(
                df_temp["jd"]-JD0, df_temp["flux_model"], marker="x",
                color=mycolor[0], label=f"LC {idx+1} model")
        ax.legend()

    # Use the same range
    if args.ylim:
        ymin, ymax = args.ylim
        for ax in fig.axes:
            ax.set_ylim([ymin, ymax])

    ax2_2.set_xlabel("Rotational Phase")
    ax2_2.set_ylabel("Relative flux")


    out = f"{args.obj}_lc{str_model}.png"
    plt.savefig(out, dpi=200)
    plt.close()
    # Plot lightcurves ========================================================
     

    # Plot phased lightcurves =================================================
    # Create 4x2 figures phased lightcurves
    if args.rotP:
        fig = plt.figure(figsize=(15, 21)) 
        ax1_1 = fig.add_axes([0.1, 0.80, 0.35, 0.18])
        ax1_2 = fig.add_axes([0.6, 0.80, 0.35, 0.18])
        ax2_1 = fig.add_axes([0.1, 0.55, 0.35, 0.18])
        ax2_2 = fig.add_axes([0.6, 0.55, 0.35, 0.18])
        ax3_1 = fig.add_axes([0.1, 0.30, 0.35, 0.18])
        ax3_2 = fig.add_axes([0.6, 0.30, 0.35, 0.18])
        ax4_1 = fig.add_axes([0.1, 0.05, 0.35, 0.18])
        ax4_2 = fig.add_axes([0.6, 0.05, 0.35, 0.18])

        if args.hour:
            rotP_sec = args.rotP*3600.
        else:
            # rotP in sec
            rotP_sec = args.rotP
        # rotP in day
        rotP_day = rotP_sec/24./3600.

        df["phase"] = [((x-JD0)/rotP_day)%1 for x in df["jd"]]
         
        for idx, ax in enumerate(
            [ax1_1, ax1_2, ax2_1, ax2_2, ax3_1, ax3_2, ax4_1, ax4_2]):
            ax.set_xlabel("Rotational Phase")
            ax.set_ylabel("Relative flux")
            df_temp = df[df["n_lc"]==idx+1]
            ax.scatter(
                df_temp["phase"], df_temp["flux"], 
                color=mycolor[1], label=f"LC {idx+1}")

            if args.lc_model:
                ax.scatter(
                    df_temp["phase"], df_temp["flux_model"], marker="x",
                    color=mycolor[0], label=f"LC {idx+1} model")
            ax.set_xlim([0.0, 1.0])
            # Ignore outliers
            #ax.set_ylim([0.5, 1.5])
            ax.legend()

        ax2_2.set_xlabel("Rotational Phase")
        ax2_2.set_ylabel("Relative flux")
        # Use the same range
        if args.ylim:
            ymin, ymax = args.ylim
            for ax in fig.axes:
                ax.set_ylim([ymin, ymax])

        # Plot lightcurves
        out = f"{args.obj}_plc{str_model}.png"
        plt.savefig(out, dpi=200)
        plt.close()
    # Plot phased lightcurves =================================================
