import os
import ftplib, urllib2
import shelve
from warnings import warn
import subprocess as sbp

import numpy as np
import pylab as pl
from pyhdf.SD import SD, SDC

from njord import nasa

class L2(object):
    """
    Read and manipulate L2 files.

    This class is designed to handle L2 LAC files from the Goddard
    Space Center, but can be used for L2 files in general.

    """

    
    def __init__(self, filename=None, jd=None, **kwargs):
        self.filename = filename
        self.jd = jd
        self.setup_grid()
        #super(L2, self).__init__(**kwargs)
        #self.add_vc()
                
    def setup_grid(self):
        """Create matrices with latitudes and longitudes for t-pos"""
        self.datadir = "data/"
        self.defaultfilename = "A2013018172500.L2_LAC_OC"
        if self.filename is not None:
            self.filename = os.path.join(self.datadir,
                                         self.filename.strip('.bz2'))
        elif self.jd is not None:
            pass
        else:
            self.filename = os.path.join(self.datadir, self.defaultfilename)
        self._try_to_unzip()
        sd = SD(self.filename, SDC.READ)
        self.llat = sd.select('latitude')[:]
        self.llon = sd.select('longitude')[:]
        self.llon[self.llon == -999.] = np.nan
        self.llat[self.llat == -999.] = np.nan
        self.j1 = 0; self.i1 = 0
        self.j2,self.i2 = self.llat.shape
        self._try_to_zip()

    def load(self,fieldname='chlor_a', nan=np.nan):
        self._l2read(fieldname, nan)

    def _l2read(self, fieldname, nan=np.nan):
        """ Read a L2 LAC file and add field to current instance"""
        if os.path.isfile(self.filename + '.npz'):
            field,base,intercept,slope = self._read_npz(fieldname)
            self.zipped = False
        else:
            self._try_to_unzip()
            field,base,intercept,slope,nanval = self._read_hdf(fieldname)
        mask = field == nanval
        if base != -999:
            self.__dict__[fieldname] = base**((slope * field) + intercept)
        else:
            self.__dict__[fieldname] = ((slope * field) + intercept)
        self.__dict__[fieldname][mask] = nan
        self._try_to_zip()
        self._try_to_npz(field, mask)

    def _read_hdf(self, fieldname='l3m_data'):
        sd = SD(self.filename, SDC.READ)
        ds = sd.select(fieldname)
        attrdict = ds.attributes()
        for d in attrdict.keys(): attrdict[d.lower()] = attrdict.pop(d)
        field      = ds[self.j1:self.j2, self.i1:self.i2].copy()
        intercept  = ds.attributes()['intercept']
        slope      = ds.attributes()['slope']
        try:
            nanval = ds.attributes()['fill']
        except:
            nanval = ds.attributes()['bad_value_scaled']
        try:
            base   = ds.attributes()['base']
        except KeyError:
            base   = -999
        dstr = sd.attributes()['Start Time'] 
        self.jd = pl.datestr2num(dstr[:4] + "-1-1") + float(dstr[4:7]) - 1
        return field,base,intercept,slope,nanval

    def _try_to_unzip(self):
        """Unzip file if exists and and is a valid bzip2 file"""
        zipfname = self.filename + '.bz2'
        if not os.path.isfile(self.filename):
            err = sbp.call(["pbzip2", "-d", zipfname])
            if err == 1:
                raise IOError, "Decompression of " + zipfname + " failed."
            self.zipped = True
        else:
            self.zipped = False

    def _try_to_zip(self):
        if self.zipped is True:
            err = sbp.call(["pbzip2", self.filename])
            if err ==1 :
                raise IOerror( "Compression of " + self.filename + " failed.")

    def _try_to_npz(self, field, mask):
        pass
        #if ((not os.path.isfile(filename + '.npz')) &
        #    ((self.llat.shape) == (self.jmt,self.imt))):
        #    self.add_ij()
        #    self._l3write_npz(filename, field[~mask],
        #                      self.imat[~mask] ,self.jmat[~mask],
        #                      base, intercept, slope)   


    def add_ij(self):
        self.imat,self.jmat = np.meshgrid(np.arange(self.i2-self.i1),
                                          np.arange(self.j2-self.j1))
        self.kdijvec = np.vstack((np.ravel(self.imat),
                                  np.ravel(self.jmat))).T



