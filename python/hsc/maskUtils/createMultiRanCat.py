#!/usr/bin/env python
#
# Jean Coupon (jean.coupon@unige.ch)
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/) and the HSC software team
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#


import numpy as np
import errno
import os

from argparse import ArgumentError


import lsst.pex.config     as pexConfig
from lsst.pipe.tasks.coaddBase import CoaddBaseTask
import lsst.afw.table as afwTable


__all__ = ["CreateMultiRanCatTask"]

class CreateMultiRanCatConfig(CoaddBaseTask.ConfigClass):

    filters   = pexConfig.Field("Name of filters to combine [default HSC-G^HSC-R^HSC-I^HSC-Z^HSC-Y]", str, "HSC-G^HSC-R^HSC-I^HSC-Z^HSC-Y")

    fileOutName = pexConfig.Field("Name of output file", str, "")
    dirOutName  = pexConfig.Field("Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH/MultiRanCat-FILTER-TRACT-PATCH.fits)", str, "")

    dustSgpFileName = pexConfig.Field("Name of output file", str, "/Users/coupon/data/SchlegelDust/SFD_dust_4096_sgp.fits")
    dustNgpFileName = pexConfig.Field("Name of output file", str, "/Users/coupon/data/SchlegelDust/SFD_dust_4096_ngp.fits")

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

class CreateMultiRanCatTask(CoaddBaseTask):
    """A Task to merge catalogs
    """

    _DefaultName = 'CreateMultiRanCat'
    ConfigClass = CreateMultiRanCatConfig

    class dustMap(object):
        """Dust map info
        """
        def __init__(self):
            pass

    def __init__(self, schema=None, *args, **kwargs):
        CoaddBaseTask.__init__(self,  *args, **kwargs)

        # ---------------------------------------------------------- #
        # for Galactic extinction until https://hsc-jira.astro.princeton.edu/jira/browse/HSC-1350 is fixed
        # ---------------------------------------------------------- #
        from   astropy.io        import ascii,fits
        import astropy.wcs       as wcs

        sFile = fits.open(self.config.dustSgpFileName)
        nFile = fits.open(self.config.dustNgpFileName)

        self.dustMap.sMap  = sFile[0].data
        self.dustMap.nMap  = nFile[0].data

        self.dustMap.sWcs = wcs.WCS(sFile[0].header)
        self.dustMap.nWcs = wcs.WCS(nFile[0].header)
        # ---------------------------------------------------------- #
        #
        # ---------------------------------------------------------- #


    def iround(self, x):
        """iround(number) -> integer
        Round a number to the nearest integer.
        From https://www.daniweb.com/software-development/python/threads/299459/round-to-nearest-integer-
        """
        return int(round(x) - .5) + (x > 0)

    def readCatalog(self, dataRef, filterName):
        """Read input catalog

        We read the input dataset provided by the 'inputDataset'
        class variable.
        """

        dataRef.dataId["filter"] = filterName

        dirIntName = dataRef.getButler().mapper.root+"/"+self.config.coaddName+"Coadd-results"
        fileInName = "{0}/{1}/{2}/{3}/ran-{1}-{2}-{3}.fits".format(dirIntName,dataRef.dataId["filter"],dataRef.dataId["tract"],dataRef.dataId["patch"])

        catalog = afwTable.SourceCatalog.readFits(fileInName, 0)
 #dataRef.get("deepCoadd_forced_src", immediate=True)
        self.log.info("Read %d sources for filter %s: %s" % (len(catalog), filterName, dataRef.dataId))
        return filterName, catalog

    def readCoadd(self, dataRef, filterName):
        """Read input coadd
        """
        dataRef.dataId["filter"] = filterName
        coadd      = dataRef.get("deepCoadd_calexp")
        self.log.info("Read coadd for filter %s: %s" % (filterName, dataRef.dataId))
        return filterName, coadd


    def getDustCorrection(self, dustMap, ra, dec):

        from astropy import units as u
        from astropy.coordinates import SkyCoord
        import astropy.wcs       as wcs

        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')

        if coord.galactic.b.degree > 0.0:
            x, y = wcs.utils.skycoord_to_pixel(coord, dustMap.nWcs,  origin=0)
            return float(dustMap.nMap[self.iround(y), self.iround(x)])
        else:
            x, y = wcs.utils.skycoord_to_pixel(coord, dustMap.sWcs,  origin=0)
            return float(dustMap.sMap[self.iround(y), self.iround(x)])


    def run(self, dataRef, selectDataList=[]):

        self.log.info("Processing %s" % (dataRef.dataId))

        filters   = self.config.filters.split("^")
        catalogs = dict(self.readCatalog(dataRef, f) for f in filters)
        coadds   = dict(self.readCoadd(dataRef, f) for f in filters)

        wcs         = coadds[filters[0]].getWcs()
        pixel_scale = wcs.pixelScale().asDegrees()*3600.0


        ref = catalogs["HSC-I"]

        # create new table table
        mergedSchema = afwTable.Schema()

        fields=[]
        # define table fields
        fields.append(mergedSchema.addField("id",               type="L", doc="Unique id"))
        fields.append(mergedSchema.addField("ra",               type="F", doc="ra [deg]"))
        fields.append(mergedSchema.addField("dec",              type="F", doc="dec [deg]"))
        fields.append(mergedSchema.addField("countInputs",      type="I", doc="Number of input single exposures for the reference filter"))
        fields.append(mergedSchema.addField("PSFDetRadius",     type="F", doc="Determinant radius for the PSF at the object position = sigma if gaussian [arcsec]"))
        fields.append(mergedSchema.addField("EB_V",             type="F", doc="Milky Way dust E(B-V) [mag]"))

        fields.append(mergedSchema.addField("isDuplicated",     type="I", doc="1 if outside the inner tract or patch"))
        fields.append(mergedSchema.addField("isOffImage",       type="I", doc="1 if in NO_DATA area or on the CCD edge"))
        fields.append(mergedSchema.addField("hasBadPhotometry", type="I", doc="1 if interpolated, saturated, suspect, has CR at center or near bright object"))
        fields.append(mergedSchema.addField("isClean",          type="I", doc="1 if none of other flags is set"))


        # create table object
        merged = afwTable.BaseCatalog(mergedSchema)

        N = len(ref)
        for i in range(N):
        #for i in range(100,110):

            # create new record
            record = merged.addNew()
            coord = ref[i].get('coord')

            # record if any of the filter is flagged as bad photometry
            for f in filters:
                hasBadPhotometry = (catalogs[f][i].get('flags.pixel.interpolated.center')) \
                                |  (catalogs[f][i].get('flags.pixel.saturated.center')) \
                                |  (catalogs[f][i].get('flags.pixel.suspect.center'))  \
                                |  (catalogs[f][i].get('flags.pixel.cr.center')) \
                                |  (catalogs[f][i].get('flags.pixel.bright.object.center'))
                if hasBadPhotometry:
                    break

            isDuplicated = not ref[i].get('detect.is-primary')
            isOffImage   = (ref[i].get('flags.pixel.offimage')) | (ref[i].get('flags.pixel.edge'))
            isClean = (not hasBadPhotometry) & (not isDuplicated) & (not isOffImage)


            # record common info from reference filter
            record.set(mergedSchema['id'].asKey(),               ref[i].get('id'))
            record.set(mergedSchema['ra'].asKey(),               coord.toFk5().getRa().asDegrees())
            record.set(mergedSchema['dec'].asKey(),              coord.toFk5().getDec().asDegrees())
            record.set(mergedSchema['countInputs'].asKey(),      ref[i].get('countInputs'))
            record.set(mergedSchema['PSFDetRadius'].asKey(),     ref[i].get("shape.sdss.psf").getDeterminantRadius()*pixel_scale)

            record.set(mergedSchema['isDuplicated'].asKey(),     int(isDuplicated))
            record.set(mergedSchema['isOffImage'].asKey(),       int(isOffImage))
            record.set(mergedSchema['hasBadPhotometry'].asKey(), int(hasBadPhotometry))
            record.set(mergedSchema['isClean'].asKey(),          int(isClean))

            # dust correction
            EB_V = self.getDustCorrection(self.dustMap, record.get(mergedSchema['ra'].asKey()), record.get(mergedSchema['dec'].asKey()))
            record.set(mergedSchema['EB_V'].asKey(), EB_V)



        # write catalog
        if self.config.fileOutName == "":

            # get output dir
            # TO DO: create new PAF
            # see /Users/coupon/local/source/hscPipe/install/DarwinX86/solvetansip/6.5.1p_hsc/python/hsc/meas/tansip/utils.py
            # and /Users/coupon/local/source/hscPipe/install/DarwinX86/pex_policy/HSC-4.0.0/tests/Policy_1.py

            if self.config.dirOutName == "" :
                dirOutName = dataRef.getButler().mapper.root+"/"+self.config.coaddName+"Coadd-results"
                self.log.info("WARNING: the output file will be written in {0:s}.".format(dirOutName))
            else:
                dirOutName = self.config.dirOutName

            fileOutName = "{0}/{1}/{2}/{3}/MultiRanCat-{1}-{2}-{3}.fits".format(dirOutName,"merged",dataRef.dataId["tract"],dataRef.dataId["patch"])


        else:
            fileOutName = self.config.fileOutName

        self.log.info("Writing {0:s}".format(fileOutName))

        self.mkdir_p(os.path.dirname(fileOutName))
        merged.writeFits(fileOutName)




        # write catalog
        #self.log.info("Writing {0:s}".format(self.config.fileOutName))
        #if self.config.fileOutName == "":
        #    fileOutName = "{0}/{1}/{2}/{3}/MultiRanCat-{2}-{3}.fits".format(self.config.dirOutName,"merged",dataRef.dataId["tract"],dataRef.dataId["patch"])
        #    self.mkdir_p(os.path.dirname(fileOutName))
        #else:
        #    fileOutName = self.config.fileOutName
        #merged.writeFits(fileOutName)

        return


    # Don't forget to overload these
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None


    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

if __name__ == '__main__':
    CreateMultiRanCatTask.parseAndRun()
