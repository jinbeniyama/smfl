#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot phase function (not so-called phase curve) used in DAMIT.

References
----------
Kaasalainen et al. 2001b, Icarus, 153, 37.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  


def get_args():
    parser = ap(description="Plot the phase function in DAMIT.")
    return parser.parse_args()
 

def main(args=None):
    if args is None:
        args = get_args()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])

    ax.set_xlabel("Phase angle [deg]")
    ax.set_ylabel("Relative brightness")

    # DAMIT phase function
    #
    # f(alpha) = 1 + A0 * exp(-alpha / D) + k * alpha
    #
    # IMPORTANT:
    #   alpha : radian
    #   D     : radian
    #   k     : 1/radian
    #
    # Plot axis is shown in degree.

    # Default values from DAMIT document
    # Example additional parameter set:
    #   Higher opposition surge and steeper phase darkening
    
    # Default values from DAMIT document
    # Additional example:
    #   Softer opposition surge and weaker phase darkening
    
    A0_list = [0.5, 0.2]
    D_list  = [0.1, 0.2]
    k_list  = [-0.5, -0.2]

    # Display in degree
    alpha_deg = np.arange(0.0, 100.0, 0.1)

    # Compute in radian
    alpha = np.radians(alpha_deg)

    for idx, (A0, D, k) in enumerate(zip(A0_list, D_list, k_list)):

        f = 1.0 + A0 * np.exp(-alpha / D) + k * alpha

        if idx == 0:
            label = (
                f"DAMIT Default "
                f"(A0, D, k)=({A0}, {D}, {k})"
            )
        else:
            label = f"(A0, D, k)=({A0}, {D}, {k})"

        ax.plot(alpha_deg, f, label=label)

    ax.legend()
    y0, y1 = ax.get_ylim()
    ax.set_ylim([0, y1])

    out = "phasefunc_in_damit.png"
    plt.savefig(out, dpi=200)


if __name__ == "__main__":
    main()
