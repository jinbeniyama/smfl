#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert shape model (.obj) created with Metasequoia
for input of convexinv and lcgenerator.
"""
from argparse import ArgumentParser as ap


if __name__ == "__main__":
    parser = ap(description="Convert shape model from Metasequoia.")
    parser.add_argument(
        "model", type=str,
        help="Original model")
    parser.add_argument(
        "out", type=str, 
        help="output filename")
    args = parser.parse_args()
   
    with open(args.model, "r") as f:
        # Skip 4 headers   
        lines = f.readlines()[4:]

        for l in lines:
            if l[0]=="#":
                # Extract
                # "# N vertices"
                # "# N elements"
                ls = l.split(" ")
                if len(ls)==3:
                    if ls[2] == "vertices\n":
                        N_v = int(ls[1])
                    elif ls[2] == "elements\n":
                        N_f = int(ls[1])
                else:
                    pass

        # Save x, y, z of vertexes
        x_list, y_list, z_list = [], [], []
        # Save combinations of facets
        f1_list, f2_list, f3_list = [], [], []
        for l in lines:
            ls = l.split()

            if len(ls)==0:
                pass
            elif (len(ls)==4) & (ls[0]=="v"):
                _, x, y, z = l.split()
                x_list.append(float(x))
                y_list.append(float(y))
                z_list.append(float(z))
            elif (len(ls)==4) & (ls[0]=="f"):
                _, f1, f2, f3 = l.split()
                f1 = f1.split("/")[0]
                f2 = f2.split("/")[0]
                f3 = f3.split("/")[0]
                f1_list.append(int(f1))
                f2_list.append(int(f2))
                f3_list.append(int(f3))
            else:
                pass

    with open(args.out, "w") as f:
        # Number of vertexes and facets
        f.write(f"{N_v} {N_f}\n")

        # x, y, z of vertexes
        for x, y, z in zip(x_list, y_list, z_list):
            f.write(f"{x} {y} {z}\n")
 
        # Combinations of facets
        for f1, f2, f3 in zip(f1_list, f2_list, f3_list):
            f.write(f"{f1} {f2} {f3}\n")
