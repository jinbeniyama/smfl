#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot shape model of asteroid while fitting with ellipsoid.

x, y, z: vertex

Reference
---------
https://notebook.community/tylerjereddy/pycon-2016/computational_geometry_tutorial
"""

import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
import matplotlib
from scipy.spatial import ConvexHull

from smfl import mycolor, fit_ellipsoid, draw_ellipsoid


if __name__ == "__main__":
    parser = ap(description="Plot shape model")
    parser.add_argument(
        "model", help="output of convexinv minkowski")
    parser.add_argument(
        "obj", type=str, help="object name")
    parser.add_argument(
        "--align", action="store_true", 
        help="Fit data to ellipsoid aligned with xyz axes")
    args = parser.parse_args()


    x, y, z = [], [], []
    with open(args.model, "r") as f:
      f = f.readlines()[1:]
      for line in f:
        line = line.split(" ")
        line = [x.strip("\n") for x in line if x!="" and x!="\n"]
        assert len(line)==3, f"Invalid input.: n_row={len(line)}"
        if float(line[0]).is_integer():
          break
        x.append(float(line[0]))
        y.append(float(line[1]))
        z.append(float(line[2]))

    
    
    fig = plt.figure(figsize=(24, 8))
    ax1 = fig.add_axes([0.1, 0.1, 0.25, 0.85], projection="3d")
    ax2 = fig.add_axes([0.4, 0.1, 0.25, 0.85], projection="3d")
    ax3 = fig.add_axes([0.7, 0.1, 0.25, 0.85], projection="3d")
    # Set plot directions with DAMIT
    # Elev from z-axis !!!
    ax1.view_init(elev=0, azim=0)
    ax2.view_init(elev=0, azim=90)
    ax3.view_init(elev=90, azim=270)

    # elev=0, azim=0
    ax1.set_title(f"{args.obj} (from +x-axis)")
    # elev=0, azim=90
    ax2.set_title(f"{args.obj} (from +y-axis)")
    # elev=90, azim=0
    ax3.set_title(f"{args.obj} (from +z-axis)")

    pts = np.array([[x,y,z] for (x,y,z) in zip(x,y,z)])
    hull = ConvexHull(pts)
    # hull_faces[i] consists of 3x3 elements
    hull_facets = hull.points[hull.simplices]
    
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    flu_triangles = Poly3DCollection(hull_facets, alpha = 1, color="gray")
    ax1.add_collection3d(flu_triangles)
    flu_triangles = Poly3DCollection(hull_facets, alpha = 1, color="gray")
    ax2.add_collection3d(flu_triangles)
    flu_triangles = Poly3DCollection(hull_facets, alpha = 1, color="gray")
    ax3.add_collection3d(flu_triangles)

    # # Plot facet
    # for x in hull_facets:
    #   X, Y, Z = x
    #   assert Fala
    #   ax1.plot_surface(X, Y, Z, cmap='jet')

 
    # Plot edge
    for (idx,s) in enumerate(hull.simplices):
      if idx==0:
        label = "convex"
      else: 
        label = None
      # Here we cycle back to the first coordinate
      s = np.append(s, s[0])  
      # ax1.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-",  color="gray", label=label)
      # ax2.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-",  color="gray", label=label)
      # ax3.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-",  color="gray", label=label)
    

    # Vertex
    for ax in [ax1, ax2, ax3]:
      ax.scatter(x, y, z, s=0.1, alpha=0.5,color="black", label=None)
      ax.set_xlabel("x")
      ax.set_ylabel("y")
      ax.set_zlabel("z")
      ax.set_box_aspect((1,1,1))
      #ax.set_zlim([-2,2])
      ax.legend()

    # Set plot range
    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()
    zmin, zmax = ax1.get_zlim()
    plot_min, plot_max = np.min([xmin, ymin, zmin]), np.max([xmax, ymax, zmax])
    for ax in [ax1, ax2, ax3]:
      ax.set_xlim([plot_min, plot_max])
      ax.set_ylim([plot_min, plot_max])
      ax.set_zlim([plot_min, plot_max])
    print(f"{plot_min}--{plot_max}")

    # Plot ellipsoid
    x, y, z = np.array(x), np.array(y), np.array(z)
    if args.align:
        # Fit with ellipsoids aligned with x, y, z-axes
        center, evecs, radii, v = fit_ellipsoid(x, y, z, True)
        str_align = "_aligned"
    else:
        # Fit with ellipsoids not aligned
        center, evecs, radii, v = fit_ellipsoid(x, y, z, False)
        str_align = ""

    a0, b0, c = radii
    if a0 > b0:
        a, b = a0, b0
    else:
        a, b = b0, a0
    x0, y0, z0 = center
    print(f"center       : x0, y0, z0 = {x0:.2f}, {y0:.2f}, {z0:.2f}")
    print(f"radii        :  a,  b,  c = {a:.2f}, {b:.2f}, {c:.2f}")

    print(f"axial ratio  :  b/a,  c/a  = {b/a:.2f}, {c/a:.2f}")

    for ax in [ax1, ax2, ax3]:
        draw_ellipsoid([0, 0, 0], radii, evecs, ax=ax, plot_axes=True, cage_color='g')

     
    out = f"{args.obj}_convexshape{str_align}.png"
    plt.savefig(out)
