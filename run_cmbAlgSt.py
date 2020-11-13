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
    n1=5000
    n2=5200
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
    piaSRT=fh['/NS/SRT/pathAtten'][n1:n2,:]
    reliabFlag=fh['NS/SRT/reliabFlag'][n1:n2,:]
    binBB=fh['NS/CSF/binBBPeak'][n1:n2,:]
    binBBT=fh['NS/CSF/binBBTop'][n1:n2,:]
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
import time
t1=time.time()
dmT=cmb.tablep2.dmj[:289]
rT=cmb.tablep2.rj[:289]
dnL=[]
dnCoeff=array([-0.00570364,  0.13319214])
dn=1.5*(dnCoeff[0]*cmb.tablep2.zkusj[:289]+dnCoeff[1])
zAd=cmb.tablep2.zkusj[:289]+10*dn
rAd=rT*10**dn

for i in range(50):
    r=i+1
    dm=(r/1.37)**(1/4.2)
    ifind = cmb.bisection2(dmT,dm)
    dn1=log10(r/rT[ifind-1])
    dnL.append([log(dm),dn1,cmb.tablep2.zkusj[ifind-1]-10*dn1])
    
for (i1,j1) in zip(a[0],a[1]):
    if j1>=12 and j1<37:
        zKu1=zKu[i1,j1,:]
        zKa1=zKa[i1,j1-12,:]
        if bzd[i1,j1]<bcf[i1,j1] and (reliabFlag[i1,j1]==1 or reliabFlag[i1,j1]==2):
            it+=1
            imu=3
            dr=0.125
            eps=1.0
            dn1d,dm1d,rrate1d,\
                zkuc,zkasim,epst,piaKu,piaKa =\
                    cmb.iter_profcv(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],\
                                  bsfc[i1,j1],zKu1,zKa1,dr,eps,imu)
            r1L.append([dprsfcRate[i1,j1],rrate1d[bcf[i1,j1]],dn1d[bcf[i1,j1]],dm1d[bcf[i1,j1]],epst,\
                        piaKu,piaSRT[i1,j1]])

print("dt=",time.time()-t1)
plt.figure()
plt.subplot(211)
plt.pcolormesh(zKu[:,24,::-1].T,cmap='jet',vmax=45)
plt.ylim(0,60)
plt.subplot(212)
plt.plot(pType[:,24])

rrate1dL=[]
for i1 in range(200):
    rrate1d=zeros(176)
    j1=24
    zKu1=zKu[i1,j1,:]
    zKa1=zKa[i1,j1-12,:]
    if pType[i1,24]==2 or pType[i1,24]==3:
        dn1d,dm1d,rrate1d,\
            zkuc,zkasim,epst,piaKu,piaKa =\
                cmb.iter_profcv(bst[i1,24],bzd[i1,24],bcf[i1,24],\
                                bsfc[i1,24],zKu1,zKa1,dr,eps,imu)
    if pType[i1,24]==1:
        if binBB[i1,j1]<=1:
            dn1d,dm1d,rrate1d,\
                zkuc,zkasim,epst,piaKu,piaKa =\
                    cmb.iter_profcv(bst[i1,24],bzd[i1,24],bcf[i1,24],\
                                bsfc[i1,24],zKu1,zKa1,dr,eps,imu)
        else:
            bbb=binBB[i1,j1]+2
            bbt=binBBT[i1,j1]
            bb=binBB[i1,j1]
            dn1d,dm1d,rrate1d,\
                zkuc,zkasim,epst,piaKu,piaKa =\
                    cmb.iter_profst(bst[i1,24],bzd[i1,24],bb,bbt,bbb,bcf[i1,24],\
                                    bsfc[i1,24],zKu1,zKa1,dr,eps,imu)
        
    rrate1dL.append(rrate1d)


rrate2d=array(rrate1dL)
rrate2d=ma.array(rrate2d,mask=rrate2d<0.001)
plt.figure()
plt.pcolormesh(rrate2d[:,::-1].T,cmap='jet')
plt.ylim(0,60)
plt.xlim(75,175)
plt.colorbar()
plt.show()
