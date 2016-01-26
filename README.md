# maskUtils

This module adds a few mask utilities to the HSC pipeline environment.

To install it, download the latest version to path/to/maskUtils
```
git clone https://github.com/jcoupon/maskUtils.git
```
and run
```
setup -v -j -r path/to/maskUtils
```

## drawRandoms.py

Creates a random catalogue for a given patch:
```
 drawRandoms.py DATADIR --rerun RERUNDIR --id tract=TRACT patch=PATCH filter=FILTER --config Nden=100  dirOutName=DIR  
```

### Options:

```Nden```: controls the number density of random points per sq. arcmin  (default is 100).

```N```: controls the number of random points per patch (default is 100000). IMPORTANT: Supersedes Nden if set. Use with caution, patches (even within a single tract) may not have the same size, so that using N=fixed would result in non-constant random point density.

```dirOutName```: Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH/ran-FILTER-TRACT-PATCH.fits (default is current directory).

```fileOutName```: the name (with full path) of the output fits file (supersedes dirOutName).

### Output:

```adjust_density```: a random number between 0 and 1 to adjust the sky density. A cut like ```adjust_density < 0.5``` will decrease the sky density by a factor of two.

```sky_mean```: Mean of sky value in 2" diameter apertures (```merge.footprint.sky == True``` objects). Constant within patch.
```sky_std``` : standard deviation of sky value in 2" diameter apertures (```merge.footprint.sky == True``` objects). Constant within patch. To compute the limiting magnitude at 5 sigma: lim_mag = -2.5 log10(5*sky_std) + 27.0

```pix_variance``` : Pixel variance at random point position.

Below is the current masking information generated by the pipeline (v 3.6.1):

| Label |	Meaning (Note that there is no order to these)|
| ------------- |-------------| -----|
| BAD	| Pixel is physically bad (a known camera defect) |
| CR |	Cosmic Ray hit  |
| CROSSTALK	| Pixel location affected by crosstalk (and corrected)  |
| EDGE   | Near the CCD edge  |
| INTERPOLATED	| Pixel contains a value based on interpolation from neighbours.  |
| INTRP| 	(same as INTERPOLATED)  |
| SATURATED	| Pixel flux exceeded full-well  |
| SAT	| (same as SATURATED)  |
| SUSPECT	| Pixel is nearly saturated. It may not be well corrected for non-linearity.  |
| UNMASKEDNAN	| A NaN occurred in this pixel in ISR (instrument signature removal - bias,flat,etc) |
| DETECTED	| Pixel is part of a source footprint (a detected source) |
| DETECTED\_NEGATIVE	 | Pixel is part of a negative source footprint (in difference image)|
| CLIPPED| 	(Coadd only) Co-addition process clipped 1 or 2 (but not more) input pixels|
| NO_DATA	| (Coadd only) Pixel has no input data (between CCDs, beyond edge of frame)|

See documentation for masks [here](http://hsca.ipmu.jp/hscsphinx/pipeline_tools.html#masks)
