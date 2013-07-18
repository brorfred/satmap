"""Generate historical data from L3 data for a given region """

import tables as td
import numpy as np
import pylab as pl
from scipy.stats import nanmean

from njord import nasa, oscar
import projmap
import figpref

miv = np.ma.masked_invalid


class L3(nasa.MODIS):

    #def __init__(self, lon1=-71,lon2=-20,lat1=35,lat2=50,res="4km",**kwargs):
    def __init__(self, lon1=-30,lon2=-15,lat1=40,lat2=50,res="4km",**kwargs):
        self.projname = "nasa.MODIS"
        super(L3, self).__init__(lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2,
                                 res=res, map_region='dimensions', **kwargs)
        jd1 = pl.datestr2num('2003-01-01')
        jd2 = pl.datestr2num('2013-05-30')+1
        self.jdvec = np.arange(jd1, jd2)

    def create_h5(self, fldname='chl'):

        h5f = h5f = td.openFile('dimensionsL3mat.h5', 'a')
        fatom = td.FloatCol()
        filtr = td.Filters(complevel=5, complib='zlib')
        crc = h5f.createCArray

        shp = (len(self.jdvec), self.llat.shape[0], self.llat.shape[1])
        mat = crc(h5f.root, fldname, fatom, shp, filters=filtr)
        for n,jd in enumerate(self.jdvec):
            self.load(fldname, jd=jd)
            mat[n,:,:] = self.chl
            print n
        llat = crc(h5f.root, 'llat', fatom,  shp[1:], filters=filtr)
        llat[:] = self.llat
        llon = crc(h5f.root, 'llon', fatom,  shp[1:], filters=filtr)
        llon[:] = self.llon
        jdvec = crc(h5f.root, 'jdvec', fatom, (shp[0],), filters=filtr)
        jdvec[:] = self.jdvec
        h5f.close()


    def stat(self, fldname='chl',scale='lin'):
        if not hasattr(self,'h5f'): self.h5open()
        h5fld = self.h5f.root.chl
        sum = np.zeros((h5fld.shape[1:]))
        cnt = np.zeros((h5fld.shape[1:]))
        var = np.zeros((h5fld.shape[1:]))
        for m in h5fld:
            mask = ~np.isnan(m)
            if scale is 'lin':
                sum[mask] += m[mask]
            else:
                sum[mask] += np.log(m[mask])
            cnt[mask] += m[mask]*0+1
        avg = sum/cnt
        for m in h5fld:
            mask = ~np.isnan(m)
            var[mask] += (m[mask]-avg[mask])**2
        std = np.sqrt(var/cnt)
        setattr(self,'%ssum%s' % (fldname,scale), sum)
        setattr(self,'%scnt%s' % (fldname,scale), cnt)
        setattr(self,'%savg%s' % (fldname,scale), avg)
        setattr(self,'%sstd%s' % (fldname,scale), std)

    def h5open(self, filename='dimensionsL3mat.h5'):
        self.h5f = td.openFile('dimensionsL3mat.h5', 'a')

    def diff(self):
        diff = zeros(chl.shape,dtype=float32)
        for t in range(chl.shape[0]-1):
            diff[t,...] = chl[t+1,...] - chl[t,...]
            print t

    def dayplot(self, jd=733773.0, djd=5):
        
        if not hasattr(self,'os'):
            self.os = oscar.Oscar(lat2=50,lat1=35,lon1=-71,lon2=-25)
        if not hasattr(self,'h5f'):
            self.h5open()
        chl = self.h5f.root.chl
        tpos = np.nonzero(self.jdvec==jd)[0][0]

        pl.clf()
        pl.subplots_adjust(hspace=0)
        pl.subplot(2,1,1)
        x,y = self.mp(self.llon,self.llat)
        self.mp.pcolormesh(x,y,miv(nanmean(np.log(chl[tpos:tpos+djd,:,:]),axis=0)))
        self.mp.nice()

        pl.subplot(2,1,2)
        self.os.load(jd=jd)
        x,y = self.mp(self.os.llon,self.os.llat)
        self.mp.contourf(x, y, miv(np.sqrt(self.os.u**2+self.os.v**2)),
                         np.arange(0,1.5,0.05))
        self.mp.nice()


    def time(self, i=500,j=250):

        if not hasattr(self,'h5f'): self.h5open()
        mat = self.h5f.root.chl[:,j-5:j+5,i-5:i+5]

        figpref.presentation()
        jd =  pl.datestr2num('2012-01-01')
        jd2 = pl.datestr2num('2013-04-30')+1
        pl.gca().xaxis.axis_date()
        pl.scatter(self.jdvec[:,np.newaxis,np.newaxis]+mat*0,mat, 5,'g')
        pl.xlim(jd,jd2)
        pl.scatter(self.jdvec,nanmean(nanmean(mat,axis=1),axis=1),20,'y')
        #setp(gca(),yscale='log')
        pl.ylim(0.01,5)
        return pl.gca()

    def all_timeseries(self):

        mp = projmap.Projmap('dimensions')
        if not hasattr(self,'h5f'): self.h5open()
        for j in np.arange(5,self.h5f.root.chl.shape[1],50):
            for i in np.arange(5,self.h5f.root.chl.shape[2],50):
                pl.clf()
                ax = self.time(i=i,j=j)
                ax.set_ylim(0.01,5)
                pl.setp(ax, yscale='log')
                bb = ax.get_position()
                ax2 = pl.axes(list(bb.p0) + [0.3,0.15])
                x,y = mp(self.llon[j-5:j+5,i-5:i+5], self.llat[j-5:j+5,i-5:i+5])
                mp.scatter(x,y,5,'c')
                mp.nice(latlabels=False, lonlabels=False)
                pl.savefig('time_log_black_%03i_%03i.pdf' % (i,j))
                print i,j


def timecube(region=None, lat1=None, lat2=None, lon1=None, lon2=None):

    if region == "BATS":
        lat1 = 30
        lat2 = 34
        lon1 = -66
        lon2 = -62
    if region == "experiment":
        lat1 = 43
        lat2 = 47
        lon1 = -26
        lon2 = -22
    ns = nasa.MODIS(lat1=lat1,lat2=lat2,lon1=lon1,lon2=lon2)
           
    jd1 = pl.datestr2num('2003-01-01')
    jd2 = pl.datestr2num('2013-05-31')
    ns.jdvec = np.arange(jd1,jd2+1)
    ns.timecube = np.zeros((len(ns.jdvec),) + ns.llat.shape)

    for n,jd in enumerate(ns.jdvec):
       ns.timecube[n,:,:] =  ns.get_field('chl',jd=jd)
       print ns.jdvec.max() -jd
    return ns
 
    
def plot_timeseries(ns,nsexp):
    pl.clf()
    pl.subplot(2,1,1)
    pl.scatter(ns.jdvec[:,np.newaxis,np.newaxis] + 
               ns.timecube[:,20:30,20:30]*0,ns.timecube[:,20:30,20:30],2,'r')
    pl.scatter(ns.jdvec,nanmean(nanmean(ns.timecube[:,20:30,20:30],axis=1),axis=1),5,'y')
    pl.legend(('Satellite Observations','Daily means'))
    pl.gca().xaxis.axis_date()
    pl.xlim(pl.datestr2num('2003-01-01'), pl.datestr2num('2013-05-31'))
    pl.setp(pl.gca(), yscale="log")
    pl.ylim(0.01,5)
    pl.title('BATS')
    pl.ylabel(r'Chl (mg m$^{-3}$ d$^{-1}$)')

    pl.subplot(2,1,2)
    pl.scatter(nsexp.jdvec[:,np.newaxis,np.newaxis] + 
               nsexp.timecube[:,20:30,20:30]*0,
               nsexp.timecube[:,20:30,20:30],2,'r')
    pl.scatter(nsexp.jdvec,
               nanmean(nanmean(nsexp.timecube[:,20:30,20:30],axis=1),axis=1),5,'y')
    pl.gca().xaxis.axis_date()
    pl.xlim(pl.datestr2num('2003-01-01'), pl.datestr2num('2013-05-31'))
    pl.setp(pl.gca(), yscale="log")
    pl.ylabel(r'Chl (mg m$^{-3}$ d$^{-1}$)')
    pl.ylim(0.01,5)
    pl.title(r'Experiment Site (45$\degree$N 24$\degree$W)')
