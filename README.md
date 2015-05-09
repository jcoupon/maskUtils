# maskUtils

This module adds a few mask utilities

For example, to create a random catalogue for a given patch:

> setup -v -j -r path/to/maskUtils
> drawRandoms.py DATADIR --rerun RERUNDIR --id tract=XXX patch=X,Y filter=HSC-X --config N=100000

N controls the number of random points per patch (default is 100000).

At the moment it uses pyfits and returns a fits table "ran.fits" in the working directory, but it would be better if it can be integrated in the pipeline product structure. Perhaps in "deepCoadd-results/HSC-X/TRACT/X,Y/ran-HSC-X-TRACT-X,Y.fits" ?
