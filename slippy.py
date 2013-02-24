"""
Create Slippymaps from numpy arrays


Reference:
    http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
"""
import os

import numpy as np
import pylab as pl
from scipy.spatial import cKDTree

class TileArray(object):

    def __init__(self, llon, llat, lon1=None, lon2=None,
                 lat1=None, lat2=None, zoom=4):
        alist = ['zoom', 'numtiles', 'llon', 'llat']
        vlist = [zoom,   2**zoom,    llon,   llat]
        for a,v in zip(alist,vlist): self.__dict__[a] = v
        self.lonvec = np.ravel(self.llon)
        self.latvec = np.ravel(self.llat)
        
        self.lat1 = self.llat.min()
        if lat1 is not None: self.lat1 = lat1
        self.lat2 = self.llat.max()
        if lat2 is not None: self.lat2 = lat2
        self.lon1 = self.llon.min()
        if lon1 is not None: self.lon1 = lon1
        self.lon2 = self.llon.max()
        if lon2 is not None: self.lon2 = lon2

        self.create_ijmats()
        self.create_llmats()

        print len(self.llon.flat)
        print len(self.clon.flat)


    def __setattr__(self, name, value):
        """Oveload __setitem__ to make sure that zoom level can change."""
        self.__dict__[name] = value
        if name == "zoom":
            self.set_zoom(value)

    def set_zoom(self, zoom):
        self.__dict__['zoom'] = zoom
        self.numtiles = 2**zoom
        self.create_ijmats()
        self.create_llmats()


    def create_ijmats(self):
        """Create arrays defining tiles at a given zoom level"""
        self.xtile1,self.ytile2 = self.deg2num(self.lon1, self.lat1)
        self.xtile2,self.ytile1 = self.deg2num(self.lon2, self.lat2)
        self.xnumtiles = self.xtile2 - self.xtile1 + 1
        self.ynumtiles = self.ytile2 - self.ytile1 + 1
        self.xtilevec = np.linspace(self.xtile1, self.xtile2+1,
                                    256*self.xnumtiles, endpoint=False)
        self.ytilevec = np.linspace(self.ytile1, self.ytile2+1,
                                    256*self.ynumtiles, endpoint=False)
        self.xtilemat,self.ytilemat = np.meshgrid(self.xtilevec,self.ytilevec)
        self.imat,self.jmat = np.meshgrid(np.arange(256*self.xnumtiles),
                                          np.arange(256*self.ynumtiles))
        self.itilemat = self.xtilemat.astype(np.int)
        self.jtilemat = self.ytilemat.astype(np.int)

    def create_llmats(self):
        """Create arrays with lon-lats for tile cells"""
        dx = self.xtilevec[1] - self.xtilevec[0]
        dy = self.ytilevec[1] - self.ytilevec[0]
        self.nwlon,self.nwlat = self.num2deg(self.xtilemat, self.ytilemat)
        self.selon,self.selat = self.num2deg(self.xtilemat+dx,
                                             self.ytilemat+dy)
        self.clon = (self.nwlon + self.selon)/2
        self.clat = (self.nwlat + self.selat)/2

    def add_kd(self, coord="tile",mask=None):
        """Generate a KD-tree objects and for the current njord instance"""
        if coord == "grid":
            if mask is None: mask = self.lonvec == self.lonvec
            xvec = self.kdivec = self.lonvec[np.ravel(mask)]
            yvec = self.kdjvec = self.latvec[np.ravel(mask)]
        else:
            xvec = self.kdlonvec = np.ravel(self.clon)
            yvec = self.kdlatvec = np.ravel(self.clat)
        self.kd = cKDTree(list(np.vstack((xvec, yvec)).T))
        #self.kdijvec = np.vstack((np.ravel(self.imat),
        #                          np.ravel(self.jmat))).T

    def reproject(self, fld, zoom=None):
        """Reproject a lat-lon vector to i-j grid coordinates"""
        if zoom is not None: self.zoom = zoom
        mask = np.ravel(~np.isnan(fld))

        self.add_kd('grid', mask=mask)
        dist,ij = self.kd.query(self.zip(np.ravel(self.clon),
                                         np.ravel(self.clat)), 1)
        tiles = self.clat * np.nan
        tiles.flat = fld.flat[mask][ij]
        tiles.flat[dist>.04] = np.nan

        self.add_kd('tile')
        dist,ij = self.kd.query(self.zip(self.lonvec[mask],
                                         self.latvec[mask]), 1)
        tilesums = np.bincount(ij, np.ravel(fld)[mask])
        tilecnts = np.bincount(ij)
        mask = tilecnts > 0
        tiles.flat[mask] = (tilesums.astype(np.float)/tilecnts)[mask]
        tiles[self.selat<self.lat1] = np.nan
        tiles[self.nwlat>self.lat2] = np.nan
        tiles[self.selon<self.lon1] = np.nan
        tiles[self.nwlon>self.lon2] = np.nan
        return tiles

    def save(self, fld, name):
        server_dir = "/Users/bror/brorfred.org/"
        pardir = "%s/tiles/%s" % (server_dir, name)
        self.safemakedirs(pardir)
        cmin = -2.5  #np.nanmin(np.log(fld))
        cmax =  1.0  #np.nanmax(np.log(fld))
        for zm in np.arange(0,8):
            print "Generating zoom level %i" % zm
            tiles = self.reproject(fld, zm)
            jsh,ish = tiles.shape
            for i in np.arange(0,ish,256):
                for j in np.arange(0,jsh,256):
                    tiledir= ("%s/%i/%i" % (pardir, self.zoom,
                                            self.itilemat[j,i]))
                    self.safemakedirs(tiledir)
                    filename = "%s/%i.png" % (tiledir, self.jtilemat[j,i])
                    pl.imsave(filename, np.log(tiles[j:j+255, i:i+255]),
                              vmin=cmin, vmax=cmax)

    def ijinterp(self,ivec,jvec, field, mask=None, nei=3, dpos=None):
        dist,ij = self.kd.query(list(np.vstack((ivec,jvec)).T), nei)
        sumvec = np.nansum(field[self.kdijvec[ij][:,:,1],
                                 self.kdijvec[ij][:,:,0]]*0+1,axis=1)
        weights = (1-dist/sumvec[:,np.newaxis])
        weights = weights/weights.sum(axis=1)[:,np.newaxis]
        fldvec = weights[:,0] * 0
        for n in np.arange(nei):
            fldvec = fldvec + field[self.kdijvec[ij[:,n]][:,1],
                                    self.kdijvec[ij[:,n]][:,0]] * weights[:,n]
        if dpos is not None:
            ipos = self.kdijvec[ij[dpos,:], 0]
            jpos = self.kdijvec[ij[dpos,:], 1]
        return fldvec

    def num2deg(self, xtile, ytile):
        """Convert from tile ID to lat-lon """
        n = self.numtiles
        lon = xtile / n * 360.0 - 180.0
        lat_rad = np.arctan(np.sinh(np.pi * (1 - 2 * ytile / n)))
        lat = np.degrees(lat_rad)
        return (lon, lat)

    def deg2num(self, lon, lat):
        """Convert from lat-lon to tile ID"""
        lat_rad = np.radians(lat)
        n = self.numtiles
        xtile = int((lon + 180.0) / 360.0 * n)
        ytile = int((1.0 - np.log(np.tan(lat_rad) + (1/np.cos(lat_rad))) /
                     np.pi) / 2.0 * n)
        return (xtile, ytile)

    def safemakedirs(self,path):
        try:
            os.makedirs(path)
        except OSError:
            pass


    def zip(self, xarr, yarr):
        """
        Transform data in a ziplike fashion. NDarrays will be flattened.

        Example:

        >>> xarr = np.array([1,2,3])
        >>> yarr = np.array([5,6,7])
        >>> self.zip(xarr, yarr)
        array([[1, 5],
               [2, 6],
               [3, 7]])

        """
        return np.vstack((np.ravel(xarr), np.ravel(yarr))).T
