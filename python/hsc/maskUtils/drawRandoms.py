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

# from lsst.obs.hsc import HscMapper

__all__ = ["DrawRandomsTask"]

class DrawRandomsConfig(CoaddBaseTask.ConfigClass):

    N           = pexConfig.Field("Number of random points per patch (supersedes Nden)", int, -1)
    Nden        = pexConfig.Field("Random number density per sq arcmin", float, 100)
    clobber     = pexConfig.Field("To overwrite existing file [default: True]", bool, True)
    #fileOutName = pexConfig.Field("Name of output file", str, "ran.fits")
    dirOutName  = pexConfig.Field("Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH/ran-FILTER-TRACT-PATCH.fits)", str, "")
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

    def makeIdFactory(self, dataRef):
        """Return an IdFactory for setting the detection identifiers

        The actual parameters used in the IdFactory are provided by
        the butler (through the provided data reference.
        """
        datasetName="MergedCoaddId"

        expBits = dataRef.get(self.config.coaddName + datasetName + "_bits")
        expId = long(dataRef.get(self.config.coaddName + datasetName))

        return afwTable.IdFactory.makeSource(expId, 64 - expBits)

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

        # print measureSourcesConfig.algorithms["flags.pixel"]


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
        table  = afwTable.SourceTable.make(self.schema, self.makeIdFactory(dataRef))
        catalog = afwTable.SourceCatalog(table)
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

            size_psf = 1.0
            try:
                shape_sdss_psf_val = psf.computeShape(afwGeom.Point2D(xy))
            except :
                pass
            else:
                record.set(shape_sdss_psf, shape_sdss_psf_val)
                size_psf = shape_sdss_psf_val.getDeterminantRadius()

            foot = afwDetection.Footprint(afwGeom.Point2I(xy), size_psf)
            record.setFootprint(foot)

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
            # try:
            #    ms.applyWithPeak(record, coadd) # need a footprint
            # except:
            #    pass

            ms.apply(record, coadd, afwGeom.Point2D(xy))

        self.setPrimaryFlags.run(catalog, skyInfo.skyMap, skyInfo.tractInfo, skyInfo.patchInfo, includeDeblend=False)

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

            fileOutName = "{0}/{1}/{2}/{3}/ran-{1}-{2}-{3}.fits".format(dirOutName,dataRef.dataId["filter"],dataRef.dataId["tract"],dataRef.dataId["patch"])


        else:
            fileOutName = self.config.fileOutName

        self.mkdir_p(os.path.dirname(fileOutName))
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

#    def check_bit(self, bitmask, bit):
#        return ((bitmask&(1<<bit))!=0)


#    def int_robust(self, s):
#        try:
#            return int(s)
#        except ValueError:
#            return 0
