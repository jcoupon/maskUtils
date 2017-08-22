#!/usr/bin/env python
#
# Jean Coupon (jean.coupon@unige.ch)
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

import healpy

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.image as afwImage
from   lsst.afw.geom import Box2D
import lsst.afw.table as afwTable
import lsst.meas.base as measBase

import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection

from lsst.pipe.tasks.coaddBase import CoaddBaseTask
from lsst.pipe.tasks.setPrimaryFlags import SetPrimaryFlagsTask

__all__ = ["DrawRandomsTask"]

class DrawRandomsConfig(CoaddBaseTask.ConfigClass):
    """configuration for drawRandoms
    """

    N = pexConfig.Field("Number of random points per patch (supersedes Nden)", int, -1)
    Nden = pexConfig.Field("Random number density per sq arcmin", float, 100)
    clobber = pexConfig.Field("To overwrite existing file [default: True]", bool, True)
    dirOutName = pexConfig.Field("Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH/ran-FILTER-TRACT-PATCH.fits)", str, "")
    fileOutName = pexConfig.Field("Name of output file (supersedes dirOutName)", str, "")
    test = pexConfig.Field("To write a test table", bool, False)
    seed = pexConfig.Field("Seed for random generator (default: based on patch id)", int, -1)
    setPrimaryFlags = pexConfig.ConfigurableField(target=SetPrimaryFlagsTask, doc="Set flags for primary source in tract/patch")
    depthMapFileName = pexConfig.Field("Name of healpix file that records full depth", str, "")

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

class DrawRandomsTask(CoaddBaseTask):
    _DefaultName = "drawRandoms"
    ConfigClass  = DrawRandomsConfig

    class depthMap(object):
        """depthMap map info in
        healpix format
        """

        def __init__(self):
            pass

    def __init__(self, schema=None, *args, **kwargs):
        CoaddBaseTask.__init__(self,  *args, **kwargs)

        if schema is None:
            schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema = schema

        self.makeSubtask("setPrimaryFlags", schema=self.schema)

        self.depthMap.map = None

        if self.config.depthMapFileName != "":
            """healpix map
            """

            self.depthMap.nest = True
            self.depthMap.map, h = healpy.fitsfunc.read_map(self.config.depthMapFileName, nest=self.depthMap.nest, h=True)
            self.depthMap.header = dict(h)
            self.depthMap.nside = self.depthMap.header['NSIDE']

    def makeIdFactory(self, dataRef):
        """ Return an IdFactory for setting the detection identifiers

        The actual parameters used in the IdFactory are provided by
        the butler (through the provided data reference).
        """

        # see /data1a/ana/hscPipe5/Linux64/afw/5.2-hsc/tests/testSourceTable.py

        datasetName="CoaddId"
        expBits = dataRef.get(self.config.coaddName + datasetName + "_bits")
        expId = long(dataRef.get(self.config.coaddName + datasetName))

        return afwTable.IdFactory.makeSource(expId, 64 - expBits)

    def run(self, dataRef, selectDataList=[]):
        """Draw randoms for a given patch
        """

        # first test if the forced-src file exists
        # do not process if the patch doesn't exist
        try:
            dataRef.get(self.config.coaddName + "Coadd_forced_src")
        except:
            self.log.info("No forced_src file found for %s. Skipping..." % (dataRef.dataId))
            return

        # verbose
        self.log.info("Processing %s" % (dataRef.dataId))

        # create a seed that depends on patch id
        # so it is consistent among filters
        if self.config.seed == -1:
            p = [int(d) for d in dataRef.dataId["patch"].split(",") ]
            numpy.random.seed(seed=dataRef.dataId["tract"]*10000+p[0]*10+ p[1])
        else:
            numpy.random.seed(seed=self.config.seed)

        # compute sky mean and sky std_dev for this patch
        # in 2" diameter apertures (~12 pixels x 0.17"/pixel)
        # import source list for getting sky objects
        sources = dataRef.get(self.config.coaddName + "Coadd_meas")
        if True:
            sky_apertures = sources['base_CircularApertureFlux_12_0_flux'][sources['merge_peak_sky']]
            select = numpy.isfinite(sky_apertures)
            sky_mean = numpy.mean(sky_apertures[select])
            sky_std  = numpy.std(sky_apertures[select])
            # NOTE: to get 5-sigma limiting magnitudes:
            # print -2.5*numpy.log10(5.0*sky_std/coadd.getCalib().getFluxMag0()[0])
        else:
            sky_mean = 0.0
            sky_std  = 0.0

        # get coadd, coadd info and coadd psf object
        coadd = dataRef.get(self.config.coaddName + "Coadd_calexp")
        psf = coadd.getPsf()
        var = coadd.getMaskedImage().getVariance().getArray()
        skyInfo = self.getSkyInfo(dataRef)

        # wcs and reference point (wrt tract)
        # See http://hsca.ipmu.jp/hscsphinx_test/scripts/print_coord.html
        # for coordinate routines.
        wcs = coadd.getWcs()
        xy0 = coadd.getXY0()

        # dimension in pixels
        dim = coadd.getDimensions()

        # define measurement algorithms
        # mostly copied from /data1a/ana/hscPipe5/Linux64/meas_base/5.3-hsc/tests/testInputCount.py
        measureSourcesConfig = measBase.SingleFrameMeasurementConfig()
        measureSourcesConfig.plugins.names = ['base_PixelFlags', 'base_PeakCentroid', 'base_InputCount', 'base_SdssShape']
        # measureSourcesConfig.plugins.names = ['base_PixelFlags', 'base_PeakCentroid', 'base_InputCount', 'base_shapeHSM_HsmPsfMoments']
        measureSourcesConfig.slots.centroid = "base_PeakCentroid"
        measureSourcesConfig.slots.psfFlux = None
        measureSourcesConfig.slots.apFlux = None
        measureSourcesConfig.slots.modelFlux = None
        measureSourcesConfig.slots.instFlux = None
        measureSourcesConfig.slots.calibFlux = None
        measureSourcesConfig.slots.shape =  None

        # it seems it is still necessary to manually add the
        # bright-star mask flag by hand
        measureSourcesConfig.plugins['base_PixelFlags'].masksFpCenter.append("BRIGHT_OBJECT")
        measureSourcesConfig.plugins['base_PixelFlags'].masksFpAnywhere.append("BRIGHT_OBJECT")

        measureSourcesConfig.validate()

        # add PSF shape
        # sdssShape_psf = self.schema.addField("shape_sdss_psf", type="MomentsD", doc="PSF xx from SDSS algorithm", units="pixel")
        # shape_sdss_psf = self.schema.addField("shape_sdss_psf", type="MomentsD", doc="PSF yy from SDSS algorithm", units="pixel")
        # shape_sdss_psf = self.schema.addField("shape_sdss_psf", type="MomentsD", doc="PSF xy from SDSS algorithm", units="pixel")

        # additional columns

        # random number to adjust sky density
        adjust_density = self.schema.addField("adjust_density", type=float, doc="Random number between [0:1] to adjust sky density", units='')

        # sky mean and variance for the entire patch
        sky_mean_key = self.schema.addField("sky_mean", type=float, doc="Mean of sky value in 2\" diamter apertures", units='count')
        sky_std_key  = self.schema.addField("sky_std", type=float, doc="Standard deviation of sky value in 2\" diamter apertures", units='count')

        # pixel variance at random point position
        pix_variance = self.schema.addField("pix_variance", type=float, doc="Pixel variance at random point position", units="flx^2")

        # add healpix map value (if healpix map is given)
        if self.depthMap.map is not None:
            depth_key = self.schema.addField("isFullDepthColor", type="Flag", doc="True if full depth and full colors at point position", units='')

        # task and output catalog
        task = measBase.SingleFrameMeasurementTask(self.schema, config=measureSourcesConfig)
        table = afwTable.SourceTable.make(self.schema, self.makeIdFactory(dataRef))
        catalog = afwTable.SourceCatalog(table)

        if self.config.N == -1:
            # to output a constant random 
            # number density, first compute 
            # the area in degree
            pixel_area = coadd.getWcs().pixelScale().asDegrees()**2
            area = pixel_area * dim[0] * dim[1]
            N = self.iround(area*self.config.Nden*60.0*60.0)
        else:
            # fixed number if random points
            N = self.config.N

        # verbose
        self.log.info("Drawing %d random points" % (N))

        # loop over N random points
        for i in range(N):
        # for i in range(100):

            # draw one random point
            x = numpy.random.random()*(dim[0]-1)
            y = numpy.random.random()*(dim[1]-1)

            # get coordinates
            radec = wcs.pixelToSky(afwGeom.Point2D(x + xy0.getX(), y + xy0.getY()))
            xy = wcs.skyToPixel(radec)

            # new record in table
            record = catalog.addNew()
            record.setCoord(radec)

            # get PSF moments and evaluate size
            #size_psf = 1.0
            #try:
            #    shape_sdss_psf_val = psf.computeShape(afwGeom.Point2D(xy))
            #except:
            #    pass
            #else:
             #   record.set(shape_sdss_psf, shape_sdss_psf_val)
             #   size_psf = shape_sdss_psf_val.getDeterminantRadius()

            # object has no footprint
            foot = afwDetection.Footprint(afwGeom.Point2I(xy), 0.0)
            peak = foot.getPeaks().addNew()
            peak.setFx(xy[0])
            peak.setFy(xy[1])
            peak.setPeakValue(0.0)
            record.setFootprint(foot)

            # draw a number between 0 and 1 to adjust sky density
            record.set(adjust_density, numpy.random.random())

            # add sky properties
            record.set(sky_mean_key, sky_mean)
            record.set(sky_std_key, sky_std)

            # add local (pixel) variance
            record.set(pix_variance, float(var[self.iround(y), self.iround(x)]))

            # required for setPrimaryFlags
            record.set(catalog.getCentroidKey(), afwGeom.Point2D(xy))

            # add healpix map value
            if self.depthMap.map is not None:
                mapIndex = healpy.pixelfunc.ang2pix(self.depthMap.nside, numpy.pi/2.0 - radec[1].asRadians(), radec[0].asRadians(), nest=self.depthMap.nest)
                record.setFlag(depth_key, self.depthMap.map[mapIndex])

        # run measurements
        task.run(catalog, coadd)

        self.setPrimaryFlags.run(catalog, skyInfo.skyMap, skyInfo.tractInfo, skyInfo.patchInfo, includeDeblend=False)

        # write catalog
        if self.config.fileOutName == "":
            if self.config.dirOutName == "" :
                fileOutName = dataRef.get(self.config.coaddName + "Coadd_forced_src_filename")[0].replace('forced_src', 'ran')
                self.log.info("WARNING: the output file will be written in {0:s}.".format(fileOutName))
            else:
                fileOutName = "{0}/{1}/{2}/{3}/ran-{1}-{2}-{3}.fits".format(self.config.dirOutName,dataRef.dataId["filter"],dataRef.dataId["tract"],dataRef.dataId["patch"])
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
