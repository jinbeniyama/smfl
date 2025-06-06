# Shape Modeling From Lightcurves (smfl)
[developer mail](mailto:beniyama@ioa.s.u-tokyo.ac.jp)

[Document in Japanese](http://www.ioa.s.u-tokyo.ac.jp/~beniyama/pdf/DAMIT_JB.pdf)

## Overview
Do shape modeling of minor bodies from lightcurves.
All what needed in this script is a text file with 
`time (jd)` and `flux` (and `fluxerr`).
After procedures 1 to 8, shape model can be created with an arbitrary software.
This repository optimizes the public code in [DAMIT](https://astro.troja.mff.cuni.cz/projects/damit/).


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
Previously I used a Python wrapper to search a sidereal period, but I realized that just doing `cat lc | period_scan ...` is easier.
```
[single lightcurve]
cat (lc) | period_scan -v (input_ps) (out_ps)
```

3. Plot sidereal period vs. chi2.
```
[single lightcurve]
plot_sidP_chi2.py (obj) (out of 2.)

[multiple lightcurve for MC technique]
plot_sidP_chi2.py (obj) (out of 2.-1) (out of 2.-2) (out of 2.-3)
```
From this plot, we get a sidereal period of the target. There are several ways to estimate uncertainty of the sidereal period such as boot-strap method, Monte-Carlo method, and 3-sigma-like threshold (See Fatka et al. 2025, A&A for details). The determination of the sidereal period is very important; if we use a different sidereal period, we get different pole solution and shape model. Please be very catreful!

4. Make input files of convexinv.
After we get a sidereal period of the target, let's constrain pole orientation next.
Input files of convexinv (e.g., `input_ci_0_-90`) are saved in `convex_input` by default.
In this process, both sidereal period and pole orientation should be fixed as far as I know.
As a test, you can make input files specifying the size of longtitude (`L`) and latitude (`M`) as follows.
The number of combination of longitude and latitude is `LxM`.
```
[period and pole fixed, test]
make_convexinv_input.py --sidP (sidereal period in hour) --Nlam (number of longitude) --Nbeta (number of latitude) --fixP --fixpole --inpdir (directory for input files)
```

For the paper, you can use the golden spiral algorithm as follows. The number of combination of longitude and latitude is `2N+1`.
```
[period and pole fixed, golden spiral algorithm]
make_convexinv_input.py --sidP (sidereal period in hour) --N_golden (number of pole) --fixP --fixpole --inpdir (directory for input files)
```


5. Do convexinv.
Input files of convexinv (e.g., `input_ci_0_-90`) in `convex_input` are used by default.
All results and intermediate files of convexinv (e.g.,`outarea_ci_0_-90`, `outlcs_ci_0_-90`, `outpar_ci_0_-90`, `res_ci_0_-90`) are saved in `convex_result` by default.
```
[Nlam and Nbeta should be the same with 4., test]
do_convexinv.py --lc (lc created in 1.) --Nlam (number of longitude) --Nbeta (number of latitude) 

[N_golden should be the same with 4., golden spiral algorithm]
do_convexinv.py --lc (lc created in 1.) --N_golden (number of pole) 
```

6. Plot pole solution with chi2.
Image is saved in the current directory.
```
plot_polesolution.py --Nlam (number of longitude) --Nbeta (number of latitude)
```
From this plot, we constrain pole orientation. Again, there are several ways to estimate uncertainty of the pole orientation such as boot-strap method, Monte-Carlo method, and 3-sigma-like threshold (See Fatka et al. 2025, A&A for details). The determination of the pole orientation is very important.

7. Create model with fixed period and pole
```
make_convexinv_input.py --sidP  (rotation period) --lam (longitude) --beta (latitude) --fixpole --fixP
do_convexinv.py --lc (lcs) --findir . --inpdir . --final
```

8. Convert model to stl and obj
```
[model, output of the process 7.]
model2stl.py (model) --out (model.stl)
model2obj.py (model) --out (model.obj)
```

9. Plot lightcurves with model curves.
```
plot_lcs_with_model.py (object name) (preprocessed lightcurve) --rotP (rotation period) 
--lc_model  (model curve, output of the process 8.)
```

10. Derive axial ratios.
```
plot_ellipsoid.py (model) (object name)
```

## Others
- Plot phase curves used in DAMIT.
```
plot_phasecurve_damit.py
```


## Installing
Please install by pip, otherwise open paths to src and smlf directories by yourself.


## Acknowledgments
I would like to express the gratitude to the people involved in [DAMIT](https://astro.troja.mff.cuni.cz/projects/damit/).
The original paper of DAMIT is 
[Ďurech et al. (2010), DAMIT: a database of asteroid models, A&A, 513, A46](https://ui.adsabs.harvard.edu/abs/2010A%26A...513A..46D).

## References
- Kaasalainen and Torppa, 2001, Icarus, Vol. 153, 24.
- Kaasalainen, Torppa, and Muinonen, 2001, Icarus, Vol. 153, 37.
- Slivan et al. 2023, Icarus, Vol. 394, 115397. (period_scan is not used.)
  - The guideline ultimately adopted was to iterate past period and pole
    convergence only until the step change in the fit relative chi-squared decreased
    to less than about 0.5%
  - The expected symmetry of an object’s two ambiguous pole solutions
    with respect to the ‘‘photometric great circle’’ (PGC) (Magnusson et al.,
    1989) was used to perform a self-consistency check on each pair
    of derived pole solutions, as a significant departure from symmetry
    would indicate the presence of some problem with the analysis
- Fatka et al. 2025, A&A, Vol. 695, A139.
  - Sidereal period is estimated with `period_scan`
  - chi-squared values are normalized such that the global minimum equaled 1
  - Using the inverse cumulative distribution function of the chi-squared distribution, 
    3-sigma threshold for a confidence level of 0.9973 and the degrees of freedom in the model is calculated
  - During this step (constraining pole orientation), sidereal period and pole prientation were fixed
- Vokrouhlicky et al. 2017, AJ, 153, 270.
  - number of parameters are 87 for two targets
  - I guess 3 parameters for spin, 3 parameters for scattering, 9x9 parameters for shape
- Hanus et al. 2023, A&A, Vol. 679, A56.
  - 3 parameters for spin, 3 parameters for scattering, 7x7 parameters for shape
- `convexinv_doc.pdf` in the DAMIT code
  - 4 parameters for scattering (`a`, `d`, `k`, and `c`)
  - parameter `c` has usually only little effect on the solution so you can fix it at, e.g., `c = 0.1`

## Dependencies
This repository is depending on `Python`, `NumPy`, `pandas`, `SciPy`, `Astropy`, `Astroquery`, `DAMIT`.
Scripts are developed on `Python 3.9.6`, `NumPy 1.26.4`, `pandas 2.2.3`, `SciPy 1.13.1`, `Astropy 6.0.1`, `Astroquery 0.4.9.post1`, `DAMIT 0.2.1`.

## LICENCE

This software is released under the XX License.
