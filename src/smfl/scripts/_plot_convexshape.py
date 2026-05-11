#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot shape model of asteroids.

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
import datetime
import seaborn as sns
import matplotlib
import warnings
from scipy.spatial import ConvexHull

warnings.simplefilter('ignore', category=matplotlib.MatplotlibDeprecationWarning)

from myio import get_filename
from myplot import mycolor

plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10

if __name__ == "__main__":
  parser = ap(description="Plot shape model")
  parser.add_argument(
    "model", help="output of convexinv minkowski")
  parser.add_argument(
    "obj", type=str, help="object name")
  args = parser.parse_args()

  filename = get_filename(args.model)

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
      print(line)
  
  fig = plt.figure(figsize=(24, 8))
  ax1 = fig.add_axes([0.1, 0.1, 0.25, 0.85], projection="3d")
  ax2 = fig.add_axes([0.4, 0.1, 0.25, 0.85], projection="3d")
  ax3 = fig.add_axes([0.7, 0.1, 0.25, 0.85], projection="3d")
  # Set plot directions with DAMIT
  ax1.view_init(elev=0, azim=0)
  ax2.view_init(elev=0, azim=270)
  ax3.view_init(elev=90, azim=270)

  # elev=0, azim=0
  ax1.set_title(f"{args.obj} (from x-axis)")
  # elev=0, azim=90
  ax2.set_title(f"{args.obj} (from y-axis)")
  # elev=90, azim=0
  ax3.set_title(f"{args.obj} (from z-axis)")

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

   
  out = f"{args.obj}_convexshape.png"
  plt.savefig(out)
