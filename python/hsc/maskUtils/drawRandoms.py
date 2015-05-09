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
import pyfits


import lsst.pex.config as pexConfig
import lsst.pipe.base  as pipeBase
import matplotlib.pyplot as pyplot
import lsst.afw.image      as afwImage
from lsst.afw.geom import Box2D


from lsst.pipe.tasks.processCoadd import ProcessCoaddTask
from lsst.pipe.tasks.coaddBase import ExistingCoaddDataIdContainer, getSkyInfo
from lsst.pipe.tasks.setPrimaryFlags import SetPrimaryFlagsTask


#class DrawRandomsConfig(pexConfig.Config):
class DrawRandomsConfig(ProcessCoaddTask.ConfigClass):
    N = pexConfig.Field("Number of random points per patch", int, 100000)
    
class DrawRandomsTask(ProcessCoaddTask):
    _DefaultName = "drawRandoms"
    ConfigClass  = DrawRandomsConfig

    def __init__(self, **kwargs):
        ProcessCoaddTask.__init__(self, **kwargs)
  

    def drawOnePoint(self, mask_array, dim, wcs, xy0, skyInfo):
        

        x = numpy.random.random()*(dim[0]-1)
        y = numpy.random.random()*(dim[0]-1)

        value = mask_array[y, x]

        radec = wcs.pixelToSky(x + xy0.getX(), y + xy0.getY())
        xy    = wcs.skyToPixel(radec)

        # inner patch
        innerFloatBBox = Box2D(skyInfo.patchInfo.getInnerBBox())
        isPatchInner = innerFloatBBox.contains(xy)

        # inner tract
        sourceInnerTractId = skyInfo.skyMap.findTract(radec).getId()       
        isTractInner = sourceInnerTractId == skyInfo.tractInfo.getId()

        ra  = radec.getLongitude().asDegrees()  
        dec = radec.getLatitude ().asDegrees()
        
        return ra, dec, value, isPatchInner, isTractInner


    def run(self, dataRef):

        self.log.info("Processing %s" % (dataRef.dataId))

        # initialize
        skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=dataRef)
        coadd = dataRef.get(self.config.coaddName + "Coadd")

        wcs = coadd.getWcs()
        xy0 = coadd.getXY0()   

        mask = coadd.getMaskedImage().getMask()
        dim = mask.getDimensions()
        mask_array = mask.getArray()


        ra           = [0.0]*self.config.N
        dec          = [0.0]*self.config.N
        value        = [0]*self.config.N
        isPatchInner = [True]*self.config.N
        isTractInner = [True]*self.config.N


        for i in range(self.config.N):
            ra[i], dec[i], value[i], isPatchInner[i], isTractInner[i] = self.drawOnePoint(mask_array, dim, wcs, xy0, skyInfo)



        c1 = pyfits.Column(name='ra',           format='D',  array=ra)
        c2 = pyfits.Column(name='dec',          format='D',  array=dec)
        c3 = pyfits.Column(name='value',        format='J',  array=value)
        c4 = pyfits.Column(name='isPatchInner', format='L',  array=isPatchInner)
        c5 = pyfits.Column(name='isTractInner', format='L',  array=isTractInner)
        
        tbhdu = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])


        tbhdu.writeto('ran.fits', clobber=True)

        return





       





    






