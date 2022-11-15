#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert rel mag flag (0) to abs mag flag (1) in lc file for convexinv.
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
warnings.simplefilter(
  "ignore", category=matplotlib.MatplotlibDeprecationWarning)
from astroquery.jplhorizons import Horizons
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy
from astropy.time import Time

from myio import get_filename
from myplot import mycolor


if __name__ == "__main__":
  parser = ap(description="Convert mag to flux")
  parser.add_argument(
    "lc_phot", help="photometric light curve in mag")
  args = parser.parse_args()
  

  filename = "lc_abs_mag.txt"
  with open(filename, "w") as f_out:
    with open(args.lc_phot) as f:
      f = f.readlines()
      for line in f:
        line_temp = line.strip("\n")
        line_temp = line_temp.split(" ")
        line_temp = [x for x in line_temp if x != ""]
        if len(line_temp) ==2:
          # mag to flux
          line_temp[1] = "1"
          line = " ".join(line_temp)
          line = f"{line}\n"

        f_out.write(line)
