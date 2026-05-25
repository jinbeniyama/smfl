#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot shape model of asteroids from an .obj file with perfect flat arrow axes.

x, y, z: vertex
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams["text.color"] = "white"

def get_args():
    parser = ap(description="Plot shape model from OBJ file with flat arrow axes")
    parser.add_argument(
        "model", type=str, 
        help="Path to the input .obj file")
    parser.add_argument(
        "out", type=str, 
        help="Output file")
    return parser.parse_args()

def load_obj(filename):
    """Parses a Wavefront .obj file to extract vertices and faces."""
    vertices = []
    faces = []
    
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            elements = line.split()
            prefix = elements[0]
            
            if prefix == "v":
                vertices.append([float(elements[1]), float(elements[2]), float(elements[3])])
            elif prefix == "f":
                face = [int(face_part.split('/')[0]) - 1 for face_part in elements[1:]]
                faces.append(face)
                
    return np.array(vertices), faces

def main(args=None):
    
    if args == None:
        args = get_args()

    if not os.path.exists(args.model):
        print(f"Error: File not found: {args.model}")
        return

    pts, faces = load_obj(args.model)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    zmin, zmax = np.min(z), np.max(z)
    plot_min, plot_max = np.min([xmin, ymin, zmin]), np.max([xmax, ymax, zmax])

    face_colors = []
    light_dir = np.array([1.0, 1.0, 2.0]) 
    light_dir = light_dir / np.linalg.norm(light_dir)
    
    triangles = []
    for face in faces:
        poly = pts[face]
        triangles.append(poly)
        
        if len(poly) >= 3:
            v1 = poly[1] - poly[0]
            v2 = poly[2] - poly[0]
            normal = np.cross(v1, v2)
            norm = np.linalg.norm(normal)
            if norm > 0:
                normal = normal / norm
            else:
                normal = np.array([0.0, 0.0, 1.0])
        else:
            normal = np.array([0.0, 0.0, 1.0])
            
        intensity = np.dot(normal, light_dir)
        intensity = (intensity + 1.0) / 2.0 
        brightness = 0.3 + 0.6 * intensity
        face_colors.append((brightness, brightness, brightness, 1.0))
    
    fig = plt.figure(figsize=(24, 8), facecolor="black")
    
    ax1 = fig.add_axes([0.05, 0.05, 0.28, 0.85], projection="3d", facecolor="black")
    ax2 = fig.add_axes([0.38, 0.05, 0.28, 0.85], projection="3d", facecolor="black")
    ax3 = fig.add_axes([0.71, 0.05, 0.28, 0.85], projection="3d", facecolor="black")
    
    # From X-axis (YZ)
    ax1.view_init(elev=0, azim=0)      
    # From Y-axis (XZ)
    ax2.view_init(elev=0, azim=270) 
    # From Z-axis (XY)
    ax3.view_init(elev=90, azim=270)  
    
    # Setting for arrow
    arrow_config = {
        ax1: [
            ((0.88, 0.50), "y", (0.93, 0.50)),
            ((0.50, 0.88), "z", (0.50, 0.93))
        ],
        ax2: [
            ((0.12, 0.50), "x", (0.07, 0.50)), 
            ((0.50, 0.88), "z", (0.50, 0.93))
        ],
        ax3: [
            ((0.12, 0.50), "x", (0.07, 0.50)), 
            ((0.50, 0.12), "y", (0.50, 0.07))  
        ]
    }


    for i, ax in enumerate([ax1, ax2, ax3]):
        # Add 3d mesh
        mesh = Poly3DCollection(triangles, facecolors=face_colors, edgecolor="none")
        ax.add_collection3d(mesh)
        
        ax.set_axis_off()
        ax.set_box_aspect((1, 1, 1))
        
        ax.set_xlim([plot_min, plot_max])
        ax.set_ylim([plot_min, plot_max])
        ax.set_zlim([plot_min, plot_max])
        
        # Plot arrows
        for end, label, text_pos in arrow_config[ax]:
            ax.annotate("", xy=end, xytext=(0.5, 0.5), xycoords='axes fraction',
                        arrowprops=dict(arrowstyle="-|>", color="white", lw=1.5, 
                                        mutation_scale=15, patchA=None, patchB=None), zorder=-999)
            ax.text2D(text_pos[0], text_pos[1], label, transform=ax.transAxes,
                      color="white", fontsize=30, ha="center", va="center", weight="bold")
            

    print(f"Plot range: {plot_min} -- {plot_max}")
     
    # Save
    out = args.out
    plt.savefig(out, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
    print(f"Saved: {out}")

if __name__ == "__main__":
    main()
