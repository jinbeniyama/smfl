#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert model created by convexinv to stl stype.

Model format is as below.
===
        1022        2040
 -0.18400128616339409       0.11578276342964949        1.5750349547097930     
 -0.25713723998959381       0.18891871725584919        1.5750349547097930     
 -0.21933486532886989        8.9921243701614292E-002   1.5810621458783960     
 -0.16848716006522421       -1.2378014512081099E-003   1.5850324663152890     
 -0.34575339049298892       0.41257529226198042        1.5528319842956879     
 ...
           1           2           3
           1           3           4
           5           6           7
           5           7           8
           9          10          11
===

Explanations
------------
1022      : Number of vertices
2040      : Number of facets
-0.18 ... : vertex x, y, z coordinates
1 2 3     : for each facet the number of vertices and the order numbers of 
            facet vertices (anticlockwise seen from outside the body)

We calculate normal vector to the facet with vertice coordinates in the script.
"""
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == "__main__":
    parser = ap(description="Convert model to stl.")
    parser.add_argument(
        "model", type=str, 
        help="Model created by convexinv")
    args = parser.parse_args()
    
    # Output filename
    out = f"{args.model}.stl"

    with open(args.model, "r") as f_in:
        lines = f_in.readlines()
        with open(out, "w") as f_out:
            # Header
            f_out.write(f"solid {args.model}-ascii\n")

            # Number of vertices and facets, 
            line0 = lines[0].split(" ")
            line0 = [int(x.strip()) for x in line0 if x != ""]
            N_v, N_f = line0
            print(f"N_vertices, N_facets = {N_v}, {N_f}")

            # Save vertices
            v_list = []
            for l in lines[1:N_v+1]:
                l = l.split(" ")
                l = [x.replace("\n", "") for x in l]
                l = [float(x) for x in l if x != ""]
                v_list.append(l)
            assert N_v == len(v_list), "Check the code."
            
            # Save facets
            f_list = []
            for l in lines[N_v+1:N_v+1+N_f]:
                l = l.split(" ")
                l = [x.replace("\n", "") for x in l]
                l = [int(x) for x in l if x != ""]
                f_list.append(l)
            assert N_f == len(f_list), "Check the code."


            # Loop for facets
            for idx,f in enumerate(f_list):
                print(f"idx facet {idx}")
                # f contains the number of vertices A, B, C
                idx_v1, idx_v2, idx_v3 = f
                # Coordinates (x,y,z) of A, B, and C
                print(f"idx1, idx2, idx3 = {idx_v1}, {idx_v2}, {idx_v3}")
                # Note: index !!
                v1 = v_list[idx_v1-1]
                v2 = v_list[idx_v2-1]
                v3 = v_list[idx_v3-1]
                
                # Vector AB and AC
                AB = [x-y for x, y in zip(v1,v2)]
                AC = [x-y for x, y in zip(v1,v3)]
                # Cross
                v_norm = np.cross(AB, AC)
                from scipy.linalg import norm
                w_norm = norm(v_norm)
                # Normalize the vector
                v_norm = v_norm/w_norm

                # Normal vector
                f_out.write(f"    facet normal {v_norm[0]} {v_norm[1]} {v_norm[2]}\n")
                f_out.write(f"        outer loop\n")
                f_out.write(f"            vertex {v1[0]} {v1[1]} {v1[2]}\n")
                f_out.write(f"            vertex {v2[0]} {v2[1]} {v2[2]}\n")
                f_out.write(f"            vertex {v3[0]} {v3[1]} {v3[2]}\n")
                f_out.write(f"        endloop\n")
                f_out.write(f"    endfacet\n")
            f_out.write(f"endsolid\n")
