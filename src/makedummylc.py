#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make a dummy lightcurve.
"""
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap


if __name__ == "__main__":
    parser = ap(description="Make dummy lightcurve.")
    parser.add_argument(
        "jd0", type=float,
        help="Starting time in jd")
    parser.add_argument(
        "jd1", type=float,
        help="Ending time in jd")
    parser.add_argument(
        "tstep", type=float,
        help="Time step in day")
    parser.add_argument(
        "--out", default="out_dummylc.txt", 
        help="output filename")
    args = parser.parse_args()
  

    jd_list = np.arange(args.jd0, args.jd1, args.tstep)
    flux = [1.0 for d in jd_list]
    fluxerr = [0.1 for d in jd_list]
    
    df = pd.DataFrame(dict(jd=jd_list, flux=flux, fluxerr=fluxerr))
    df.to_csv(args.out, sep=" ", index=False)
