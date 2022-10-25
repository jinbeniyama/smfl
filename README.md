# Shape Modeling From Lightcurves (smfl)
[developer mail](mailto:beniyama@ioa.s.u-tokyo.ac.jp)

[Document in Japanese](http://www.ioa.s.u-tokyo.ac.jp/~beniyama/pdf/DAMIT_JB.pdf)

## Overview

Do shape modeling of minor bodied from lightcurves.
All what needed in this script is a text file with 
`time (jd)` and `flux` (and `fluxerr`).
After procedures 1 to 8, shape model can be created with an arbitrary software.



## Procedures
1. Format lightcurves with aspect data from JPL ephemerides.
```
[for periodic analysis with MC technique]
format4convexinv.py (lc1) --jpl (jpl1) --N_mc (N) --obj (obj) --out (lcs_all)

[for shape modeling]
format4convexinv.py (lc1) (lc2) (lc3) --jpl (jpl1) (jpl2) (jpl3) --obj (obj) --out (lcs_all)
```

2. Search sidereal period with period_scan.
The result of `period_scan` is saved in the current directory.
```
[single lightcurve]
search_sidP.py (lc) --inp (inputfile) --out (out)

[multiple lightcurve for MC technique]
search_sidP.py (lc1) (lc2) --inp (inputfile) --out (out)
```

3. Plot sidereal period vs. chi2.
```
[single lightcurve]
plot_sidP_chi2.py (obj) (out of 2.)

[multiple lightcurve for MC technique]
plot_sidP_chi2.py (obj) (out of 2.-1) (out of 2.-2) (out of 2.-3)
```

4. Make input files of convexinv.
Input files of convexinv (e.g., input_ci_0_-90) are saved in `convex_input` by default.
```
[pole fixed]
make_convexinv_input.py --sidP (sidereal period in hour) 
--Nlam (number of longitude) --Nbeta (number of latitude)  --fixpole
```

5. Do convexinv.
Input files of convexinv (e.g., input_ci_0_-90) in `convex_input` are used by default.
All results and intermediate files of convexinv (e.g., outarea_ci_0_-90, outlcs_ci_0_-90, outpar_ci_0_-90, res_ci_0_-90) are saved in `convex_result` by default.
```
[Nlam and Nbeta should be the same with 4.]
do_convexinv.py --Nlam (number of longitude) --Nbeta (number of latitude) 
--lc (lc created in 1.)
```

6. Plot pole solution with chi2.
Image is saved in the current directory.
```
plot_polesolution.py --Nlam (number of longitude) --Nbeta (number of latitude) --norm
```

7. Convert model to stl.
```
[in prep.]
model2stl.py (model)
```


8. Plot lightcurves with model curves.
```
[in prep.]
plot_lcs_with_model.py
```

## Installing
Please install by pip, otherwise open paths to src and smlf directories by yourself.
## Usage

## Note

## Dependencies

This library is depending on `NumPy`, `SciPy`, `SEP`, `Astropy` 
and `Astroquery`.
Scripts are developed on `Python 3.7.10`, `NumPy 1.19.2`, `SciPy 1.6.1`,
`SEP 1.0.3`, `Astropy 4.2` and `Astroquery 0.4.1`.

## LICENCE

This software is released under the XX License.
