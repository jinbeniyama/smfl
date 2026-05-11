#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert model created by convexinv to .obj file.
"""
from argparse import ArgumentParser as ap
import numpy as np


def convert2obj(model, out):
    """
    Convert model to obj file.

    Parameters
    ----------
    model : str
        model created by convexinv
    out : str
        output file name (.obj)
    """
    with open(model, 'r') as f:
        lines = f.readlines()

    header = lines[0].split()
    n_vertices, n_faces = int(header[0]), int(header[1])

    vertices = lines[1:1 + n_vertices]
    faces = lines[1 + n_vertices:1 + n_vertices + n_faces]

    with open(out, 'w') as obj:
        for v in vertices:
            x, y, z = v.strip().split()
            obj.write(f"v {x} {y} {z}\n")
        for f in faces:
            i, j, k = f.strip().split()
            #obj.write(f"f {int(i)+1} {int(j)+1} {int(k)+1}\n")
            obj.write(f"f {int(i)} {int(j)} {int(k)}\n")

if __name__ == "__main__":
    parser = ap(description="Convert model to obj.")
    parser.add_argument(
        "model", type=str, 
        help="Model created by convexinv")
    parser.add_argument(
        "--out", type=str, default="model.obj",
        help="Output model")
    args = parser.parse_args()
    
    convert2obj(args.model, args.out)
