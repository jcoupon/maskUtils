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
import pyfits as fits

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
    test        = pexConfig.Field("To write a test table", bool, False)
    
class DrawRandomsTask(ProcessCoaddTask):
    _DefaultName = "drawRandoms"
    ConfigClass  = DrawRandomsConfig

    def __init__(self, **kwargs):
        ProcessCoaddTask.__init__(self, **kwargs)

    def run(self, dataRef):

        if self.config.test:
            self.testTable(read=False) 
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

        # create arrays   
        ra           = [0.0]*self.config.N
        dec          = [0.0]*self.config.N
        isPatchInner = [True]*self.config.N
        isTractInner = [True]*self.config.N
        flags        = [[True]*self.config.N for i in range(len(keys))]

        # loop over N random points
        for i in range(self.config.N):
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

        # create table object
        tbhdu = fits.BinTableHDU.from_columns(cols)

        # write table
        tbhdu.writeto(self.config.fileOutName, clobber=True)
      
        return

    def check_bit(self, bitmask, bit):
        return ((bitmask&(1<<bit))!=0)

    def drawOnePoint(self, mask_array, dim, wcs, xy0, skyInfo):
        
        x = numpy.random.random()*(dim[0]-1)
        y = numpy.random.random()*(dim[1]-1)

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


    def testTable(self, read=False):

        # to do
        # this table cannot be read by Topcat and I have no clue why

        # create table
        schema = afwTable.Schema()

        # define table fields
        fields = [schema.addField("ra", type="F", doc="ra")]
        fields.append(schema.addField("dec", type="F", doc="dec"))
        fields.append(schema.addField("isPatchInner", type="Flag", doc="True if inside patch inner area"))
        fields.append(schema.addField("isTractInner", type="Flag", doc="True if inside tract inner area"))

        # create table object
        catalog = afwTable.BaseCatalog(schema)

        # fill the table
        record = catalog.addNew()
        record.set(fields[0], 1.0)
        record.set(fields[1], 2.0)
        for i in range(2, 4):
            record.set(fields[i], True)

        # save the table
        catalog.writeFits("test.fits")

        if read:

            #FILE = pyfits.open("/Users/coupon/data/HSC/SSP/rerun/CLAUDS/deepCoadd-results/HSC-Y/1/5,5/src-HSC-Y-1-5,5.fits")
            FILE = pyfits.open("test.fits")
            d    = FILE[1].data

            print FILE.info()
            print FILE[0].header

            FILE.close()

        return






       





    






