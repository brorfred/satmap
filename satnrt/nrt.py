import os
import ftplib, urllib2
import shelve
from warnings import warn
import subprocess as sbp
import sqlite3
import ConfigParser
import json

import numpy as np
import pylab as pl
from pyhdf.SD import SD, SDC
import requests

from njord import nasa
import l2
import yrday
import slippy

class Base(object):
    """Base class for handling Goddard/NASA NRT files"""

    def _load_presets(self, filepref, kwargs):
        """Read and parse the config file"""
        self.basedir =  os.path.dirname(os.path.abspath(__file__))
        cfg = ConfigParser.ConfigParser()
        files = ["%s/%s.cfg" % (os.curdir, filepref),
                 "%s/.%s.cfg" % (os.path.expanduser("~"), filepref),
                 "%s/%s.cfg" % (self.basedir, filepref)]
        for fnm in files:
            cfg.read(fnm)
            if self.projname in cfg.sections():
                self.config_file = fnm
                break
        else:
            raise NameError('No subscription configured')

        def splitkey(key, val):
            if key in kwargs.keys():
                self.__dict__[key] = kwargs[key]
                del kwargs[key]
            else:
                self.__dict__[key] = val
            
        for key,val in cfg.items(self.projname):
            try:
                splitkey(key, json.loads(val))
            except ValueError:
                splitkey(key, val)

class Nrt(Base):
    """
    Read and manipulate NRT files.

    This class is designed to handle L2 LAC files from the MODIS
    Near Real TIME (NRT) system. Individual files are read using
    NORD's L2 filer parsing class and then mapped to a global 4km
    grid using  nearest neighbour lookup with KDTree. Both indivi-
    dual fields and full scenes can be saved as NPZ vectors for fast
    and size efficient storage.

    The NRT service can provide two types of files:
       OC - Ocean color fields such as Chl, POC, kdPAR, and k490.
      SST - Sea Surface Temperature fields.

    """
    def __init__(self, projname="dim", verbose=True, **kwargs):
        self.projname = projname
        self.verbose  = verbose
        self._load_presets('l2nrt',kwargs)
        self.ns   = nasa.MODIS(res=self.resolution)
        self.llat = self.ns.llat.copy()
        self.llon = self.ns.llon.copy()
        self.ft = Ftp(verbose=self.verbose)
        self.ns.add_ij()

    def refresh(self, fieldname='chlor_a'):
        self.ft.refresh()
        for jd in np.unique(self.ft.jdvec[self.ft.jdvec>0].astype(np.int)):
            self.load(fieldname, jd)
            
    def load(self, fieldname="chlor_a", jd=None):
        """Load NRT files and reproject to L3 grid"""
        ft = Ftp(verbose=self.verbose)
        self.jd = self.ft.jdvec.astype(np.int).max() if jd is None else jd
        self.chl = self.ns.llat * np.nan
        npzfilename = os.path.join(self.npzdir,
                    "NRT%04i%03i_%s.npz" % (self.yr,self.yd,fieldname))
        flist = self.ft.get_latest() if jd is None else self.ft.get_byjd(jd)
        exists = self._try_to_read_npz(npzfilename)
        if exists == len(flist):
            return
        else:
            self._merge_l2_files(flist)
            mask = ~np.isnan(self.chl)
            self._write_npz(npzfilename, self.chl[mask],
                            self.ns.imat[mask], self.ns.jmat[mask], flist)
            
    def _write_npz(self, npzfilename, dvec, ivec ,jvec, flist):
        """Create an npz file with field data as vectors""" 
        kwargs = {'dvec' : dvec,
                  'ivec' : ivec.astype(np.int16),
                  'jvec' : jvec.astype(np.int16),
                  'jd'   : self.jd, 'res' : self.resolution,
                  'flist': flist}
        np.savez_compressed(npzfilename, **kwargs)

    def _try_to_read_npz(self, npzfilename):
        """Create an npz file with field data as vectors""" 
        if not os.path.exists(npzfilename): return False
        self.chl = self.llat * np.nan
        fH = np.load(npzfilename)
        self.chl[fH['jvec'],fH['ivec']] = fH['dvec']
        return len(fH['flist'])

    def _merge_l2_files(self,flist):
        for filename in flist:
            if self.verbose: print "    Merging " + filename
            lt = l2.L2(os.path.join(self.hdfdir, filename))
            lt.load(fieldname='chlor_a')
            mask = ~np.isnan(lt.llat + lt.llon + lt.chlor_a)
            if len(lt.llon[mask]) == 0: continue
            i,j = self.ns.ll2ij(lt.llon[mask], lt.llat[mask])
            self.chl[j,i] = lt.chlor_a[mask]
            self.chl[self.chl<0.01] = np.nan
        self.jd = lt.jd

    def slippy(self, jd):
        sl = slippy.TileArray(self.llon, self.llat, maxzoom=8,
                              lat1=35, lat2=55, lon1=-81, lon2=-10, 
                              cmin=np.log(0.01), cmax=np.log(10))

        dtstr = pl.num2date(jd).strftime('%Y-%m-%d')
        sl.save(self.chl, 'dimensions/%s' % dtstr, 
                '/Users/bror/brorfred.org/')


    @property
    def yr(self):
        return pl.num2date(self.jd).year
    
    @property
    def yd(self):
        return yrday.yd(np.array(self.jd)).astype(np.int) - 1 


#class Shelvesql(object):
    

class Ftp(Base):

    def __init__(self, projname="dim", verbose=True, **kwargs):
        self.projname = projname
        self.verbose  = verbose
        self._load_presets('l2nrt',kwargs)
        self.fstat = shelve.open(self.statusfile, writeback=True)

    def download(self, filename):
        """Download a file from GSFC's website"""
        if not filename in self.fstat.keys(): self.fstat[filename] = {}
        filename = os.path.basename(filename)
        url = self.fileurl % filename
        try:
            response = urllib2.urlopen(url)
        except:
            warn("File not found on the server.\n tried %s" %
                 url, UserWarning)
            return False
        if self.fstat[filename]['downloaded'] is False:
            output = open(os.path.join(self.hdfdir, filename), 'wb')
            output.write(response.read())
            output.close()
            self.fstat[filename]['downloaded'] = True
            if self.verbose: print "Downloading " + filename
            if 'L2_LAC_' in filename:
                jd = (pl.datestr2num(filename[1:5] + '-1-1') +
                      int(filename[5:8]) + float(filename[8:10])/24 +
                      float(filename[10:12])/1440) - 1
                self.fstat[filename]['jd'] = jd
            else:
                self.fstat[filename]['jd'] = -999
        self.fstat.sync()     

    def get_filelist_from_server(self):
        """Download filelist from server"""
        if self.verbose: print "Downloading filelist from server"
        params = {"subID":self.nrtid, "addurl":1, "results_as_file":1}
        req = requests.get(self.searchurl, params=params)
        self.filelist = req.text.split()
        self.newfilelist = []
        for fn in self.filelist:
            basefn = str(os.path.basename(fn))
            if not basefn in self.fstat.keys():
                self.fstat[basefn] = {'on_server':True,
                                      'staged':False,
                                      'downloaded':False,
                                      'jd':-999}
                self.newfilelist.append(basefn)
        self.fstat.sync()     

    def refresh(self):
        """Refresh list of available files and download new from server"""
        self.get_filelist_from_server()
        for filename in self.newfilelist:
            self.fstat[filename]['staged'] = True
            self.download(filename)
        if len(self.newfilelist) == 0 : print "No new files."
        self.newfilelist = []
        self.fstat.sync()     

    def get_latest(self, type='OC'):
        """Return names of most recent files"""
        fnvec = np.array(self.fstat.keys())
        jdvec = self.jdvec.astype(np.int)
        mask  = np.array([type in s for s in fnvec])
        return fnvec[mask][jdvec[mask] == jdvec[mask].max()]

    def get_byjd(self, jd, type='OC'):
        """Return names of files from a specified julian date"""
        fnvec = np.array(self.fstat.keys())
        jdvec = self.jdvec.astype(np.int)
        mask  = np.array([type in s for s in fnvec])
        return fnvec[mask][jdvec[mask] == jd]


    @property
    def jdvec(self):
        return np.array([self.fstat[d]['jd'] for d in self.fstat])
        
#success = bz2.UncompressFile("hamlet.xml.bz2","hamlet.xml")
