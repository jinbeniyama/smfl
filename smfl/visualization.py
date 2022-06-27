#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

mymark = ["o", "^", "x", "D", "+", "v", "<", ">", "h", "H"]

mycolor = [
    "#AD002D", "#1e50a2", "#69821b", "#f055f0", "#afafb0", 
    "#0095b9", "#89c3eb", "red", "blue", "orange", "pink"] 

def plotstyle(n):
    """
    Return n-th plot style.
    """
    n_col = len(mycolor)
    n_mark = len(mymark)
    n_total = n_col*n_mark
    assert n < n_total, f"Large index detected {n}."

    color = mycolor[int(n%n_col)]
    marker = mymark[int(np.floor(n/n_col))]
    return color, marker


def figure4lc(nline):
    """
    Create and return fig and axes (6x4) for lightcurve plot.
    If nline < 6, return a small sized figure.
    """
    # Absolute values
    N_raw            = 4
    N_line           = nline
    figsize_per_line = 2.5
    figsize_xlabel   = 0.05
    figsize_line     = figsize_per_line*N_line + figsize_xlabel
    figsize_raw      = 10
    print(f"N_raw, N_line = {N_raw}, {N_line}")

    # Relative values
    ## space for x label
    offset_xlabel = (figsize_xlabel/figsize_line)

    # Offset
    offset_width  = 0.0625
    offset_height = 0.03*(6/N_line)

    raw_width   = 1/N_raw*0.98
    line_height = 1/N_line*0.98

    # Figure size in each region
    fig_width  = raw_width*0.7
    fig_height = line_height*0.7

    axes = []
    n_fig = 0
    
    fig = plt.figure(figsize=(figsize_raw, figsize_line))
    for idx_line, line in enumerate(range(N_line)):
        for idx_raw, raw in enumerate(range(N_raw)):
            x0 = offset_width + idx_raw*raw_width
            y0 = offset_xlabel + offset_height + (N_line-idx_line-1)*line_height
            x1 = fig_width 
            y1 = fig_height
            ax = fig.add_axes([x0, y0, x1, y1]) 
            axes.append(ax)
            n_fig += 1
        if idx_line == (N_line-1):
            return fig, axes
    return fig, axes
