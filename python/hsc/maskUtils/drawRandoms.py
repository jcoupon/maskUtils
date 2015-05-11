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

import lsst.pex.config   as pexConfig
import lsst.pipe.base    as pipeBase
import matplotlib.pyplot as pyplot
import lsst.afw.image    as afwImage
from   lsst.afw.geom     import Box2D
import lsst.afw.table    as afwTable

from lsst.pipe.tasks.processCoadd    import ProcessCoaddTask
from lsst.pipe.tasks.coaddBase       import ExistingCoaddDataIdContainer, getSkyInfo
from lsst.pipe.tasks.setPrimaryFlags import SetPrimaryFlagsTask


#class DrawRandomsConfig(pexConfig.Config):
# inherit from ProcessCoaddTask to use patch names in command line
class DrawRandomsConfig(ProcessCoaddTask.ConfigClass):
    N           = pexConfig.Field("Number of random points per patch", int, 100000)
    fileOutName = pexConfig.Field("Name of output file", str, "ran.fits")
    test        = pexConfig.Field("To write a test table", "Flag", False)
    
class DrawRandomsTask(ProcessCoaddTask):
    _DefaultName = "drawRandoms"
    ConfigClass  = DrawRandomsConfig

    def __init__(self, **kwargs):
        ProcessCoaddTask.__init__(self, **kwargs)

    def check_bit(self, bitmask, bit):
        return ((bitmask&(1<<bit))!=0)


    def drawOnePoint(self, mask_array, dim, wcs, xy0, skyInfo):
        
        x = numpy.random.random()*(dim[0]-1)
        y = numpy.random.random()*(dim[0]-1)

        bitmask = mask_array[y, x]

        radec = wcs.pixelToSky(x + xy0.getX(), y + xy0.getY())
        xy    = wcs.skyToPixel(radec)

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


    def testTable(self):

        # create table
        schema = afwTable.Schema()

        # define table fields
        fields = [schema.addField("ra", type="F", doc="ra")]
        fields.append(schema.addField("dec", type="F", doc="dec"))
        fields.append(schema.addField("isPatchInner", type="I", doc="True if inside patch inner area"))
        fields.append(schema.addField("isTractInner", type="I", doc="True if inside tract inner area"))

        # create table object
        table = afwTable.BaseCatalog(schema)

        record = table.addNew()
        record.set(fields[0], 1.0)
        record.set(fields[1], 1.0)
        record.set(fields[2], True)
        record.set(fields[3], True)

        table.writeFits("test.fits")

        return

    def run(self, dataRef):

        if self.config.test:
            self.testTable()    
            return


        # verbose
        self.log.info("Processing %s" % (dataRef.dataId))

        # get coadd and coadd info 
        coadd   = dataRef.get(self.config.coaddName + "Coadd")
        skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRef)

        # wcs and reference point (wrt tract)
        wcs = coadd.getWcs()
        xy0 = coadd.getXY0()  

        # extract mask and convert into array
        mask = coadd.getMaskedImage().getMask()
        dim  = mask.getDimensions()
        mask_array = mask.getArray()

        # get mask labels and bit values
        mask_labels = mask.getMaskPlaneDict()
        keys        = mask_labels.keys()
        values      = mask_labels.values()

        # create table that will store random points and mask values
        schema = afwTable.Schema()

        # define table fields
        fields = [schema.addField("ra", type="F", doc="ra")]
        fields.append(schema.addField("dec", type="F", doc="dec"))
        fields.append(schema.addField("isPatchInner", type="I", doc="True if inside patch inner area"))
        fields.append(schema.addField("isTractInner", type="I", doc="True if inside tract inner area"))
        for key in keys:
            fields.append(schema.addField(key, type="I", doc=key))

        # create table object
        table = afwTable.BaseCatalog(schema)

        # loop over N random points
        for i in range(self.config.N):
            ra, dec, bitmask, isPatchInner, isTractInner = self.drawOnePoint(mask_array, dim, wcs, xy0, skyInfo)

            record = table.addNew()
            record.set(fields[0], ra)
            record.set(fields[1], dec)
            record.set(fields[2], isPatchInner)
            record.set(fields[3], isTractInner)

            # set individual mask values
            for j in range(len(values)):
                flag = self.check_bit(bitmask, values[j])
                if flag:
                    toto = 1
                else:
                    toto = 0
                record.set(fields[4+j], toto)


        table.writeFits(self.config.fileOutName)
      
        return





       





    






