# Shape Modeling From Lightcurves (smfl)
[developer mail](mailto:beniyama@ioa.s.u-tokyo.ac.jp)

## Overview

Do shape modeling of minor bodied from lightcurves.


### Procedure
All what needed in the procedure is time (jd), flux (, and fluxerr) text!

1. Format lightcurves with JPL ephemerides.
`format4convexinv.py`

2. Search sidereal period with period_scan.
`search_sidP.py`

3. Plot sidereal period vs. chi2.
`plot_sidP_chi2.py`

4. Make input files of convexinv.
`make_convexinv_input.py`

5. Conduct convexinv.
`do_convexinv.py`

6. Plot pole solution with chi2.
`plot_polesolution.py`

7. Convert model to stl.
`model2stl.py`

Then, shape model can be created with an arbitrary software.

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
