#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test of Golden spiral algorithm.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

def golden_spiral_spherical(N):
    """
    """
    phi = (1 + np.sqrt(5)) / 2
    golden_angle = 2 * np.pi / phi**2

    i = np.arange(N)
    theta = golden_angle * i
    z = 1 - 2*i/(N - 1)                
    r = np.sqrt(1 - z*z)        

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    lam = np.arctan2(y, x)          
    beta = np.arcsin(z)             

    lam_deg = np.degrees(lam)
    # From 0 to 360 deg
    lam_deg = [x%360 for x in lam_deg]
    beta_deg = np.degrees(beta)
    
    return lam_deg, beta_deg


def golden_spiral_G10(N):
    """Generate evenly distributed points on a sphere using a symmetric golden spiral.

    Parameters
    ----------
    N : int
        Half the number of points minus one. 
        The total number of output points will be 2N + 1.

    Returns
    -------
    lon : ndarray of shape (2N + 1,)
        Longitudes of the points in degrees, in the range [0, 360].

    lat : ndarray of shape (2N + 1,)
        Latitudes of the points in degrees, in the range [-90, 90].


    Reference
    ---------
    González 2010, Math Geosci (2010) 42: 49–64.
    """

    phi = (1 + np.sqrt(5)) / 2  # golden ratio
    i = np.arange(-N, N + 1)
    lat = np.arcsin(2 * i / (2 * N + 1)) * 180 / np.pi

    # fractional part of i * golden ratio, mapped to [0, 1)
    frac = np.mod(i * phi, 1.0)
    # map to [0, 360)
    lon = frac * 360  
    return lon, lat


def mean_angular_spacing(lon, lat):
    """Estimate the mean angular spacing between points on a sphere in degrees.

    Parameters
    ----------
    lon : ndarray
        Longitudes of points in degrees.
    lat : ndarray
        Latitudes of points in degrees.

    Returns
    -------
    mean_spacing : float
        Mean angular spacing (degrees) between nearest neighbor points on the sphere.
    """
    # Convert to Cartesian coordinates on the unit sphere
    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    xyz = np.vstack([x, y, z]).T

    # Use KDTree to find nearest neighbors
    tree = cKDTree(xyz)
    # first one is itself, second is nearest neighbor
    dist, _ = tree.query(xyz, k=2)  

    # Convert chord distance to angle (arc length)
    angular_dist_rad = 2 * np.arcsin(dist[:, 1] / 2)
    angular_dist_deg = np.rad2deg(angular_dist_rad)

    return np.mean(angular_dist_deg)


def plot_on_mollweide(lam_deg, beta_deg):
    """Mollweide.
    """
    plt.figure(figsize=(8, 4.5))
    ax = plt.subplot(111)
    ax.scatter(lam_deg, beta_deg, s=5, color='darkblue')

    #ax = plt.subplot(111, projection="mollweide")
    #lam_rad = np.radians(lam_deg)
    #beta_rad = np.radians(beta_deg)
    #ax.scatter(lam_rad, beta_rad, s=5, color='darkblue')

    ax.grid(True)
    ax.set_xlabel(r"$\lambda$ [rad]")
    ax.set_ylabel(r"$\beta$ [rad]")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    N = 2000

    # N points
    #lam, beta = golden_spiral_spherical(N)
    
    # 2N + 1 points
    lam, beta = golden_spiral_G10(N)
    wid = mean_angular_spacing(lam, beta)
    print(f"width = {wid:.2f}")

    plot_on_mollweide(lam, beta)

