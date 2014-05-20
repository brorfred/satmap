"""
Create Slippymaps from numpy arrays


Reference:
    http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
"""
import os
import fnmatch
import urllib
import zipfile

import numpy as np
import pylab as pl
from scipy.spatial import cKDTree

from paramiko import SSHClient, SSHConfig
from scp import SCPClient

class TileArray(object):
    """Class to generate off-line slippy tiles from zoomable maps"""
    def __init__(self, llon, llat, **kwargs):
        """Setup the class instance"""
        for key in ['maxzoom','lon1','lon2', 'lat1', 'lat2',
                    'cmin', 'cmax', 'cmap']:
            if not key in kwargs.keys():
                self.__dict__[key] = None
            else:
                self.__dict__[key] = kwargs[key]
        self.lonvec = np.ravel(llon)
        self.latvec = np.ravel(llat)

        if self.lat1 is None:
            self.lat1 = max([self.latvec.min(), -85.0511])
        else:
            self.lat1 = max([self.lat1, -85.0511])

        if self.lat2 is None:
            self.lat2 = min([self.latvec.max(), 85.0511])
        else:
            self.lat2 = min([self.lat2, 85.0511])

        if self.lon1 is None:
            self.lon1 = self.lonvec.min()
        else:
            self.lon1 = self.lon1

        if self.lon2 is None:
            self.lon2 = self.lonvec.max()
        else:
            self.lon2 = self.lon2

        self.zoom = 0

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

    def save(self, fld, name, webdir="./tiles", maxzoom=8, zip=False, scp=False):
        pardir = "%s/%s" % (webdir, name)
        self.safemakedirs(pardir)
        cmin = self.cmin if self.cmin is not None else np.nanmin(fld) 
        cmax = self.cmax if self.cmax is not None else np.nanmax(fld) 
        if self.maxzoom is not None: maxzoom = self.maxzoom
        for zm in np.arange(maxzoom+1):
            print "Generating zoom level %i" % zm
            tiles = self.reproject(fld, zm)
            jsh,ish = tiles.shape
            for i in np.arange(0,ish-1,256):
                for j in np.arange(0,jsh-1,256):
                    tiledir= ("%s/%i/%i" % (pardir, self.zoom,
                                            self.itilemat[j,i]))
                    self.safemakedirs(tiledir)
                    filename = "%s/%i.png" % (tiledir,
                                              self.jtilemat[j,i])
                    pl.imsave(filename, tiles[j:j+255, i:i+255],
                              vmin=cmin, vmax=cmax)
                    print filename
        if zip:
            self.zipdir(pardir+".zip", pardir, basedir=webdir + "/tiles/")
            if scp:
                self.scp(pardir+".zip", "dimzip")

            
    def zipdir(self,zipfilename, dir, basedir=""):

        dir = os.path.abspath(dir)
        basedir = os.path.abspath(basedir)
        print zipfilename
        def find_files(directory, pattern):
            for root, dirs, files in os.walk(directory):
                for basename in files:
                    if fnmatch.fnmatch(basename, pattern):
                        filename = os.path.join(root, basename)
                        yield filename
    
        with  zipfile.ZipFile(zipfilename,'w', zipfile.ZIP_DEFLATED) as zip:
            for fn in find_files(dir, "*"):
                zip.write(fn, fn.replace(basedir, ''))

    def scp(self, filename, remote_path):
        config = SSHConfig()
        config.parse(open(os.path.expanduser('~/.ssh/config')))
        o = config.lookup('geodata')
        ssh_client = SSHClient()
        ssh_client.load_system_host_keys()
        ssh_client.connect(o['hostname'], username=o['user'])
        scp = SCPClient(ssh_client.get_transport())
        scp.put(filename, remote_path=remote_path)

                    
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
        

    def create_static_files(self):

        py = """
        var map = L.map('map').setView([37, -46], 4);
 
        L.tileLayer('http://{s}.tile.cloudmade.com/BC9A493B41014CAABB98F0471D759707/997/256/{z}/{x}/{y}.png', {
        maxZoom: %i,
        attribution: 'Map data &copy; <a href="http://openstreetmap.org">
        OpenStreetMap</a>
        contributors, <a href="http://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>,
        Imagery (c) <a href="http://cloudmade.com">CloudMade</a>'
        }).addTo(map);

        L.tileLayer('http://{url}/tiles/{name}/{z}/{x}/{y}.png', {
        name: 'test',
        url: 'brorfred.org',
        }).addTo(map);
        """


    def bgtileurl(self, id="BC9A493B41014CAABB98F0471D759707"):
        url = "http://tile.cloudmade.com/%s/997/256/" % id
        return url
        
    def download_bgtiles(self, webdir="./", maxzoom=8):
        """Download necessary background tiles for local use"""
        pardir = "%s/tiles/%s" % (webdir, "background")
        self.safemakedirs(pardir)
        for zm in np.arange(maxzoom+1):
            print "Generating zoom level %i" % zm
            for i in np.arange(0,2**zm):
                tiledir= ("%s/%i/%i" % (pardir, zm, i))
                self.safemakedirs(tiledir)
                for j in np.arange(0,2**zm):
                    url = "%s/%i/%i/%i.png" % (self.bgtileurl(), zm,i,j)
                    filename = "%s/%i.png" % (tiledir, j)
                    urllib.urlretrieve(url, filename)
                    print url
        
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




def maketiles(name, llon, llat, field, maxzoom=8, **kwargs):
    """Create neccesary files for a slippy map
    
    Creates a hiarchy or folders that contanis 256x256 size png's according 
    to the slippy standard used by google maps. 

    Parameters
    ----------
       name : Name of project/basename for folders
       llon : Matrix of t-pos longitudes. Same shape as field.
       llat : Matrix of t-pos latitudes. Same shape as field.
      field : Field to be shown on map

    Optional Parameters
    -------------------
    maxzoom : Max zoom-level of map. 0 is the entire glob. increases by 2**x
       lon1 : Cutoff for smallest longitude
       lon2 : Cutoff for largest longitude
       lat1 : Cutoff for smallest latitude
       lat2 : Cutoff for largest latitude
       cmin : Lower level for colorbar
       cmax : Upper level for colorbar 
       cmap : Alternative colormap
    """
    for key in ['lon1', 'lon2', 'lat1', 'lat2', 'cmin', 'cmax', 'cmap']:
        if not key in kwargs.keys():
            kwargs[key] = None

    tl = TileArray(llon, llat, **kwargs)
    tl.save(field, name, maxzoom=maxzoom)
