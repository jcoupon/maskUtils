# maskUtils

This module adds a few mask utilities

For example, to create a random catalogue for a given patch:

```
 setup -v -j -r path/to/maskUtils
 drawRandoms.py DATADIR --rerun RERUNDIR --id tract=XXX patch=X,Y filter=HSC-X --config N=100000
```

N controls the number of random points per patch (default is 100000).

At the moment it uses pyfits and returns a fits table "ran.fits" in the working directory, but it would be better if it can be integrated in the pipeline product structure. Perhaps in "deepCoadd-results/HSC-X/TRACT/X,Y/ran-HSC-X-TRACT-X,Y.fits" ?

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
