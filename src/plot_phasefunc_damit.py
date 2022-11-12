#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot phase function (not so-called phase curve) used in DAMIT.

References
----------
Kaasalainen et al. 2001b, Icarus, 153, 37.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  


if __name__ == "__main__":
    parser = ap(description="Plot the phase function in DAMIT.")
    args = parser.parse_args()
 

    fig = plt.figure(figsize=(8, 6))

    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    ax.set_xlabel("Phase angle [deg]")
    ax.set_ylabel("Relative brightness")
    
    # Default values
    #   A0, D, k = 0.5, 0.1, -0.5 (DAMIT document, convexinv_doc.pdf)
    A0_list = [0.5]
    D_list = [0.1]
    k_list = [-0.5]
    alpha = np.arange(0, 50, 0.1)

    for idx, (A0,D,k) in enumerate(zip(A0_list, D_list, k_list)):
        f = A0*np.exp(-alpha/D) + k*alpha + 1
        if idx == 0:
            label = f"DAMIT Default (A0, D, k) = ({A0}, {D}, {k})"
        else:
            label = f"(A0, D, k) = ({A0}, {D}, {k})"
        ax.plot(alpha, f, label=label)
        ax.legend()

    out = f"phasefunc_in_damit.png"
    plt.savefig(out, dpi=200)
