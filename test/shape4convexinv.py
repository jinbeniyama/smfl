#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert shape model of Ryugu (SHAPE_SFM_49k_v20180804.dat) to that 
for input of convexinv and lcgenerator.
"""
from argparse import ArgumentParser as ap


if __name__ == "__main__":
    parser = ap(description="Fotmat lightcurves for convex inversion.")
    parser.add_argument(
        "model", type=str,
        help="Original model")
    parser.add_argument(
        "out", type=str, 
        help="output filename")
    args = parser.parse_args()
   
    
    with open(args.model, "r") as f:
        lines = f.readlines()
        # Number of vertexes
        N_v = int(lines[0].split()[0])
        # Number of facets
        N_f = int(lines[N_v+1].split()[0])

        # Save x, y, z of vertexes
        x_list, y_list, z_list = [], [], []
        for l in lines[1:N_v+1]:
            _, x, y, z = l.split()
            x_list.append(float(x))
            y_list.append(float(y))
            z_list.append(float(z))

        # Save combinations of facets
        f1_list, f2_list, f3_list = [], [], []
        for l in lines[N_v+2:]:
            _, f1, f2, f3 = l.split()
            f1_list.append(int(f1))
            f2_list.append(int(f2))
            f3_list.append(int(f3))
    with open(args.out, "w") as f:
        # Number of vertexes and facets
        f.write(f"{N_v} {N_f}\n")

        # x, y, z of vertexes
        for x, y, z in zip(x_list, y_list, z_list):
            f.write(f"{x} {y} {z}\n")
 
        # Combinations of facets
        for f1, f2, f3 in zip(f1_list, f2_list, f3_list):
            f.write(f"{f1} {f2} {f3}\n")
        
