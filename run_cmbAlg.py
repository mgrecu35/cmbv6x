fname='2A.GPM.DPR.V8-20180723.20180118-S125407-E142641.022105.V06A.HDF5'
fname='2A.GPM.DPR.V8-20180723.20180215-S134018-E151250.022541.V06A.HDF5'
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import combAlg as cmb
#from julia.api import Julia
#jl=Julia()#compiled_modules=False)
from numpy import *
cmb.mainfortpy()
cmb.initp2()
fh=Dataset(fname)
nx=24
iplot=0
for i in range(2,3):
    n1=400*i
    n2=400*(i+1)
    zKu=fh['NS/PRE/zFactorMeasured'][n1:n2,:,:]
    zKa=fh['MS/PRE/zFactorMeasured'][n1:n2,:,:]
    lon=fh['NS/Longitude'][n1:n2,:]
    lat=fh['NS/Latitude'][n1:n2,:]
    bzd=fh['NS/VER']['binZeroDeg'][n1:n2,:]
    bcf=fh['NS/PRE']['binClutterFreeBottom'][n1:n2,:]
    zKu=ma.array(zKu,mask=zKu<10)
    zKa=ma.array(zKa,mask=zKa<10)
    bsfc=fh['NS/PRE/binRealSurface'][n1:n2,:]
    bst=fh['NS/PRE/binStormTop'][n1:n2,:]
    pType=(fh['NS/CSF/typePrecip'][n1:n2,:]/1e7).astype(int)
    dprsfcRate=fh['NS/SLV/precipRateESurface'][n1:n2,:]
    #
    if iplot==1:                                                
        plt.figure(figsize=(12,6))
        plt.subplot(211)
        plt.title('Ku')
        plt.pcolormesh(lon[:,24],arange(176)*0.125,zKu[:,24,::-1].T,\
                       vmin=10,vmax=40,cmap="jet")
        plt.plot(lon[:,nx],(175-bzd[0:,nx])*0.125)
        plt.ylim(0,6)
        plt.colorbar()
        plt.subplot(212)
        plt.title('Ka')
        plt.pcolormesh(lon[:,24],arange(176)*0.125,zKa[:,12,::-1].T,\
                       vmin=10,vmax=40,cmap="jet")

        plt.ylim(0,6)
        plt.colorbar()

it=0
a=nonzero(pType==2)
r1L=[]
for (i1,j1) in zip(a[0],a[1]):
    if j1>=12 and j1<37:
        zKu1=zKu[i1,j1,:]
        zKa1=zKa[i1,j1-12,:]
        if bzd[i1,j1]<bcf[i1,j1]:
            it+=1
            imu=3
            dr=0.125
            eps=1.0
            dn1d,dm1d,rrate1d,\
                zkuc,zkasim,epst =\
                    cmb.iter_prof(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],\
                                  bsfc[i1,j1],zKu1,zKa1,dr,eps,imu)
            r1L.append([dprsfcRate[i1,j1],rrate1d[bcf[i1,j1]]])
if iplot==1:                                                
    plt.figure(figsize=(11,6))
    plt.subplot(211)
    nx=24
    plt.title('Ku')
    plt.pcolormesh(lon[:200,nx],arange(176)*0.125,zKu[:200,nx,::-1].T,\
                   vmin=10,vmax=40,cmap="jet")
    plt.plot(lon[:200,nx],(175-bzd[0:200,nx])*0.125)
    plt.ylabel('Height')
    plt.ylim(0,6)
    plt.colorbar()
    plt.subplot(212)
    plt.title('Ka')
    plt.pcolormesh(lon[:200,nx],arange(176)*0.125,zKa[:200,nx-12,::-1].T,\
                   vmin=10,vmax=40,cmap="jet")
    plt.xlabel('Longitude')
    plt.ylabel('Height')
    plt.ylim(0,6)
    plt.colorbar()
    plt.savefig('orbit22541_sp1.png')

    plt.figure(figsize=(11,6))
    plt.subplot(211)
    nx=24
    plt.title('Ku')
    plt.pcolormesh(lon[200:,nx],arange(176)*0.125,zKu[200:,nx,::-1].T,\
                   vmin=10,vmax=40,cmap="jet")
    plt.plot(lon[200:,nx],(175-bzd[200:,nx])*0.125)
    plt.ylabel('Height')
    plt.ylim(0,6)
    plt.colorbar()
    plt.subplot(212)
    plt.title('Ka')
    plt.pcolormesh(lon[200:,nx],arange(176)*0.125,zKa[200:,nx-12,::-1].T,\
                   vmin=10,vmax=40,cmap="jet")
    plt.xlabel('Longitude')
    plt.ylabel('Height')
    plt.ylim(0,6)
    plt.colorbar()
    plt.savefig('orbit22541_sp2.png')


plt.figure(figsize=(12,8))
for nx1 in range(153,155):
    plt.subplot(1,2,nx1-152)
    plt.plot(zKu[nx1,24,:],176*0.125-arange(176)*0.125)
    plt.plot(zKa[nx1,12,:],176*0.125-arange(176)*0.125)
    plt.ylim(0,4)
    plt.xlim(10,30)
    plt.xlabel('dBZ')
    plt.ylabel('Height(km)')
    plt.title('Scan# %3.3i'%(nx1+n1))

plt.tight_layout(pad=3)
plt.margins(0.2)
plt.savefig('zProfs22541cv.png')

plt.figure(figsize=(12,8))
for nx1 in range(131,135):
    plt.subplot(2,2,nx1-130)
    plt.plot(zKu[nx1,24,:],176*0.125-arange(176)*0.125)
    plt.plot(zKa[nx1,12,:],176*0.125-arange(176)*0.125)
    plt.ylim(0,4)
    plt.xlim(10,35)
    plt.title('Scan# %3.3i'%(nx1+n1))
    plt.xlabel('dBZ')
    plt.ylabel('Height(km)')

plt.tight_layout(pad=3)
plt.margins(0.2)
plt.savefig('zProfs22541.png')
