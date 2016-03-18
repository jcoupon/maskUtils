#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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

import numpy
import errno
import os

import lsst.pex.config   as pexConfig
import lsst.pipe.base    as pipeBase
import lsst.afw.image    as afwImage
from   lsst.afw.geom     import Box2D
import lsst.afw.table    as afwTable
import lsst.meas.algorithms as measAlg

import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection

from lsst.pipe.tasks.coaddBase import CoaddBaseTask
from lsst.pipe.tasks.setPrimaryFlags import SetPrimaryFlagsTask

__all__ = ["DrawRandomsTask"]

class DrawRandomsConfig(CoaddBaseTask.ConfigClass):

    N           = pexConfig.Field("Number of random points per patch (supersedes Nden)", int, -1)
    Nden        = pexConfig.Field("Random number density per sq arcmin", float, 100)
    clobber     = pexConfig.Field("To overwrite existing file [default: True]", bool, True)
    #fileOutName = pexConfig.Field("Name of output file", str, "ran.fits")
    dirOutName  = pexConfig.Field("Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH/ran-FILTER-TRACT-PATCH.fits)", str, ".")
    fileOutName = pexConfig.Field("Name of output file (supersedes dirOutName)", str, "")
    test        = pexConfig.Field("To write a test table", bool, False)
    seed        = pexConfig.Field("Seed for random generator (default: based on patch id)", int, -1)

    setPrimaryFlags = pexConfig.ConfigurableField(target=SetPrimaryFlagsTask, doc="Set flags for primary source in tract/patch")

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

class DrawRandomsTask(CoaddBaseTask):
    _DefaultName = "drawRandoms"
    ConfigClass  = DrawRandomsConfig

    def __init__(self, schema=None, *args, **kwargs):
        #ProcessCoaddTask.__init__(self, **kwargs)
        CoaddBaseTask.__init__(self,  *args, **kwargs)

        if schema is None:
            schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema = schema

        self.makeSubtask("setPrimaryFlags", schema=self.schema)


    def run(self, dataRef, selectDataList=[]):
        """
        Draw randoms for a given patch
        See http://hsca.ipmu.jp/hscsphinx_test/scripts/print_coord.html
        for coordinate routines.
        """

        # verbose
        self.log.info("Processing %s" % (dataRef.dataId))

        # create a seed that depends on patch id
        # so it is consistent among filters
        if self.config.seed == -1:
            p = [int(d) for d in dataRef.dataId["patch"].split(",") ]
            numpy.random.seed(seed=dataRef.dataId["tract"]*10000+p[0]*10+ p[1])
        else:
            numpy.random.seed(seed=self.config.seed)

        # for sky objects
        sources = dataRef.get(self.config.coaddName + "Coadd_meas")

        # get coadd, coadd info and coadd psf object
        #coadd = butler.get('calexp', dataRef.dataId)
        coadd = dataRef.get(self.config.coaddName + "Coadd_calexp")
        psf   = coadd.getPsf()
        var   = coadd.getMaskedImage().getVariance().getArray()

        skyInfo = self.getSkyInfo(dataRef)
        #skyInfo = self.getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRef)

        # wcs and reference point (wrt tract)
        wcs  = coadd.getWcs()
        xy0  = coadd.getXY0()

        # dimension in pixels
        dim  = coadd.getDimensions()

        # create table and catalog
        # copied from meas_algorithms/HSC-3.9.0/tests/countInputs.py
        measureSourcesConfig = measAlg.SourceMeasurementConfig()
        measureSourcesConfig.algorithms.names = ["flags.pixel", "centroid.naive", "countInputs"]

        measureSourcesConfig.algorithms["flags.pixel"].center.append("BRIGHT_OBJECT")
        measureSourcesConfig.algorithms["flags.pixel"].any.append("BRIGHT_OBJECT")

        # see /Users/coupon/local/source/hscPipe/install/DarwinX86/meas_algorithms/HSC-4.0.0/tests/variance.py for variance measurement
        # measureSourcesConfig.algorithms.names = ["flags.pixel", "centroid.naive", "countInputs", "variance"]
        # measureSourcesConfig.algorithms["variance"].mask = ["BAD", "SAT"]

        measureSourcesConfig.slots.centroid = "centroid.naive"
        measureSourcesConfig.slots.psfFlux = None
        measureSourcesConfig.slots.apFlux = None
        measureSourcesConfig.slots.modelFlux = None
        measureSourcesConfig.slots.instFlux = None
        measureSourcesConfig.slots.calibFlux = None
        measureSourcesConfig.slots.shape =  None

        measureSourcesConfig.validate()

        # add PSF shape
        shape_sdss_psf = self.schema.addField("shape_sdss_psf", type="MomentsD", doc="PSF moments from SDSS algorithm", units="Pixels")

        # add random number to adjust sky density
        adjust_density = self.schema.addField("adjust_density", type=float, doc="Random number between [0:1] to adjust sky density", units="unitless")

        # add PSF size
        # PSF_size       = self.schema.addField("PSF_size", type=numpy.float32, doc="Size of the PSF from shape_sdss_psf (=sigma of gaussian)", units="Pixels")

        # add sky mean and sky std_dev in 2" diameter apertures
        if True:
            sky_apertures = sources["flux.aperture"][sources["merge.footprint.sky"]]
            sky_mean = numpy.mean(sky_apertures[:,2])
            sky_std  = numpy.std(sky_apertures[:,2])
        else:
            sky_mean = 0.0
            sky_std  = 0.0
        sky_mean_key = self.schema.addField("sky_mean", type=float, doc="Mean of sky value in 2\" diamter apertures", units="flux")
        sky_std_key  = self.schema.addField("sky_std", type=float, doc="Standard deviation of sky value in 2\" diamter apertures", units="flux")

        pix_variance = self.schema.addField("pix_variance", type=float, doc="Pixel variance at random point position", units="flux^2")

        # to get 5-sigma limiting magnitudes:
        # print -2.5*numpy.log10(5.0*sky_std/coadd.getCalib().getFluxMag0()[0])

        ms      = measureSourcesConfig.makeMeasureSources(self.schema)
        catalog = afwTable.SourceCatalog(self.schema)
        measureSourcesConfig.slots.setupTable(catalog.getTable())

        if self.config.N == -1:
            # constant random number density
            # Compute the area in degree
            pixel_area = coadd.getWcs().pixelScale().asDegrees()**2
            area       = pixel_area * dim[0] * dim[1]
            N          = self.iround(area*self.config.Nden*60.0*60.0)
        else:
            # fixed number if random points
            N = self.config.N

        # verbose
        self.log.info("Drawing %d random points" % (N))

        # loop over N random points
        for i in range(N):

            # draw a random point
            x = numpy.random.random()*(dim[0]-1)
            y = numpy.random.random()*(dim[1]-1)

            # get coordinates
            radec = wcs.pixelToSky(afwGeom.Point2D(x + xy0.getX(), y + xy0.getY()))
            xy    = wcs.skyToPixel(radec)

            # add record in table
            record = catalog.addNew()
            record.setCoord(radec)

            # get PSF moments and evaluate size
            try:
                shape_sdss_psf_val = psf.computeShape(afwGeom.Point2D(xy))
            except :
                pass
            else:
                record.set(shape_sdss_psf, shape_sdss_psf_val)

            # size = psf.computeShape(afwGeom.Point2D(xy)).getDeterminantRadius()
            # record.set(PSF_size, size)

            # looks like defining footprint isn't necessary
            # foot = afwDetection.Footprint(afwGeom.Point2I(xy), width)
            # record.setFootprint(foot)

            # add sky properties
            record.set(sky_mean_key, sky_mean)
            record.set(sky_std_key, sky_std)

            # add local variance
            record.set(pix_variance, float(var[self.iround(y), self.iround(x)]))

            # draw a number between 0 and 1 to adjust sky density
            record.set(adjust_density, numpy.random.random())

            # required for setPrimaryFlags
            record.set(catalog.getCentroidKey(), afwGeom.Point2D(xy))

            # do measurements (flagging and countINputs)
            ms.apply(record, coadd, afwGeom.Point2D(xy))

        self.setPrimaryFlags.run(catalog, skyInfo.skyMap, skyInfo.tractInfo, skyInfo.patchInfo, includeDeblend=False)

        # write catalog
        if self.config.fileOutName == "":
            fileOutName = "{0}/{1}/{2}/{3}/ran-{1}-{2}-{3}.fits".format(self.config.dirOutName,dataRef.dataId["filter"],dataRef.dataId["tract"],dataRef.dataId["patch"])
            self.mkdir_p(os.path.dirname(fileOutName))
        else:
            fileOutName = self.config.fileOutName
        catalog.writeFits(fileOutName)

        # to do. Define output name in init (not in paf) and
        # allow parallel processing
        # write sources
        # if self.config.doWriteSources:
        #   dataRef.put(result.sources, self.dataPrefix + 'src')

        return

    # Overload these if your task inherits from CmdLineTask
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

    def iround(self, x):
        """
        iround(number) -> integer
        Round a number to the nearest integer.
        From https://www.daniweb.com/software-development/python/threads/299459/round-to-nearest-integer-
        """
        return int(round(x) - .5) + (x > 0)


    """
    *******************************************************************
    *******************************************************************
                        BELOW IS NOT USED ANYMORE
    *******************************************************************
    *******************************************************************
    """


    def run_Old(self, dataRef, selectDataList=[]):

        import pyfits as fits

        if self.config.test:
            self.testTable(read=True)
            return

        # verbose
        self.log.info("Processing %s" % (dataRef.dataId))

        # get coadd and coadd info
        coadd   = dataRef.get(self.config.coaddName + "Coadd")

        #coadd.setPsf(None)
        #print coadd.getPsf()

        skyInfo = self.getSkyInfo(dataRef)
        #skyInfo = self.getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRef)

        # wcs and reference point (wrt tract)
        wcs = coadd.getWcs()
        xy0 = coadd.getXY0()

        # extract mask and convert into array
        mask = coadd.getMaskedImage().getMask()
        dim  = mask.getDimensions()
        mask_array = mask.getArray()

        # Compute the area in degree
        pixel_area =  coadd.getWcs().pixelScale().asDegrees()**2
        area       = pixel_area * dim[0] * dim[1]

        # print "total area = ", area

        if self.config.N  == -1:
            N = self.iround(area*self.config.Nden*60.0*60.0)
        else:
            N = self.config.N

        # get mask labels and bit values
        mask_labels = mask.getMaskPlaneDict()
        keys        = mask_labels.keys()
        values      = mask_labels.values()

        # initialize arrays
        ra           = [0.0]*N
        dec          = [0.0]*N
        isPatchInner = [True]*N
        isTractInner = [True]*N
        flags        = [[True]*N for i in range(len(keys))]

        # loop over N random points
        for i in range(N):
            ra[i], dec[i], bitmask, isPatchInner[i], isTractInner[i] = self.drawOnePoint(mask_array, dim, wcs, xy0, skyInfo)

            for j, value in enumerate(values):
                flags[j][i] = self.check_bit(bitmask, value)

        # column definition
        cols = fits.ColDefs([fits.Column(name="ra",           format="E", array=ra)])
        cols.add_col(        fits.Column(name="dec",          format="E", array=dec))
        cols.add_col(        fits.Column(name="isPatchInner", format="L", array=isPatchInner))
        cols.add_col(        fits.Column(name="isTractInner", format="L", array=isTractInner))
        for j, key in enumerate(keys):
            cols.add_col(        fits.Column(name=key,            format="L", array=flags[j]))

        # create table object - method depends on pyfits version
        pyfits_v = [self.int_robust(v) for v in (fits.__version__).split('.')]
        if (pyfits_v[0] > 3) | (pyfits_v[0] == 3 & pyfits_v[1] >= 3):
            tbhdu = fits.BinTableHDU.from_columns(cols)
        else:
            tbhdu = fits.new_table(cols)

        # add info in header
        header = fits.Header()
        header['AREA'] = (area, "in sq. degrees")

        # write table
        hdulist = fits.HDUList([fits.PrimaryHDU(header=header), tbhdu])

        #tbhdu.writeto(self.config.fileOutName, clobber=self.config.clobber)
        hdulist.writeto(self.config.fileOutName, clobber=self.config.clobber)

        return

    def check_bit(self, bitmask, bit):
        return ((bitmask&(1<<bit))!=0)


    def drawOnePointOLD(self, mask_array, dim, wcs, xy0, skyInfo):

        x = numpy.random.random()*(dim[0]-1)
        y = numpy.random.random()*(dim[1]-1)

        bitmask = mask_array[y, x]

        radec = wcs.pixelToSky(x + xy0.getX(), y + xy0.getY())
        xy    = wcs.skyToPixel(radec)


        #self.exp.getWcs().pixelToSky(afwGeom.Point2D(self.x, self.y))

        # inner patch
        innerFloatBBox = Box2D(skyInfo.patchInfo.getInnerBBox())
        isPatchInner   = innerFloatBBox.contains(xy)

        # inner tract
        sourceInnerTractId = skyInfo.skyMap.findTract(radec).getId()
        isTractInner       = sourceInnerTractId == skyInfo.tractInfo.getId()

        # ra and dec in decimal degrees
        ra  = radec.getLongitude().asDegrees()
        dec = radec.getLatitude ().asDegrees()

        return ra, dec, bitmask, isPatchInner, isTractInner


    def testTable(self, read=False):

        #import lsst.afw.geom as afwGeom


        #print afwGeom.Point2D(100.23, 100)


        #exit(-1)

        # to do
        # this table cannot be read by Topcat and I have no clue why


#        schema = afwTable.SourceTable.makeMinimalSchema()
#        task = measAlg.SourceMeasurementTask(schema, config=config)
#        source = catalog.addNew()
#        source.set("id", 12345)



        # create table
#        schema = afwTable.Schema()
        schema = afwTable.SourceTable.makeMinimalSchema()

        # define table fields
        fields = [schema.addField("ra", type="F", doc="ra")]
        fields.append(schema.addField("dec", type="F", doc="dec"))
        fields.append(schema.addField("isPatchInner", type="Flag", doc="True if inside patch inner area"))
        fields.append(schema.addField("isTractInner", type="Flag", doc="True if inside tract inner area"))

        # create table object
#        catalog = afwTable.BaseCatalog(schema)
        catalog = afwTable.SourceCatalog(schema)


        # fill the table
        record = catalog.addNew()
#        help(record)


        record.set(fields[0], 1.0)
        record.set(fields[1], 2.0)
        for i in range(2, 4):
            record.set(fields[i], True)

        for source in record:
            print source.get("isPatchInner")

        # save the table
        catalog.writeFits("test.fits")

        if read:

            #FILE = fits.open("/Users/coupon/data/HSC/SSP/rerun/CLAUDS/deepCoadd-results/HSC-Y/1/5,5/src-HSC-Y-1-5,5.fits")
            FILE = fits.open("test.fits")
            d    = FILE[1].data

            print FILE.info()
            print FILE[0].header

            FILE.close()

        return


    def int_robust(self, s):
        try:
            return int(s)
        except ValueError:
            return 0
