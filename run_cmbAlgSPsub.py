fname='2A.GPM.DPR.V8-20180723.20180118-S125407-E142641.022105.V06A.HDF5'
fname='2A.GPM.DPR.V8-20180723.20180215-S134018-E151250.022541.V06A.HDF5'
fname='../ORO/DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20180416-S001234-E002112.023465.V06A.HDF5'
fname='2A-CS-Tasmania.GPM.DPR.V8-20180723.20200619-S194636-E194813.035844.V06A.HDF5'
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import combAlg as cmb
#from julia.api import Julia
#jl=Julia()#compiled_modules=False)
from numpy import *
cmb.mainfortpy()
cmb.initp2()

import glob
fs=glob.glob("../ORO/DPR-CS/2A*.Ku*HDF5")
fs=glob.glob("SP/2A*SP*HDF5")
icount=0
iplot=1
zKuL=[]
zKaL=[]
hL=[]
r1L=[]
ibbL=[]
rd1L=[]
for f in sorted(fs)[7:8]:
    f=f.replace("Ku","DPR")
    fh=Dataset(f)
    lat=fh['NS/Latitude'][:,24]
    lon=fh['NS/Longitude'][:,24]
    n1=0
    n2=900
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
    h0=fh['NS/VER/heightZeroDeg'][n1:n2,:]
    htop=fh['NS/PRE/heightStormTop'][n1:n2,:]
    epsDPR=fh['NS/SLV/epsilon'][n1:n2,:,:]
    a=nonzero(pType==2)
    b=nonzero((h0[a]-1200)*(h0[a]-4500)<0)
    icount+=len(a[0][b])
    if len(a[0])>0:
        print(icount,h0[a].min(),h0[a].max())
    else:
        continue
    c=nonzero(a[1][b]==24)
    n=zKu.shape[0]
    n3=int(n/3)
    plt.figure(figsize=(12,8))
    plt.suptitle('Orbit '+f.split('.')[-3])
    for i1 in range(1,3):
        plt.subplot(2,1,i1)
        plt.pcolormesh(arange(i1*n3,i1*n3+n3),arange(176)*0.125,\
                       zKu[i1*n3:i1*n3+n3,24,::-1].T,vmin=0,vmax=45,cmap="jet")
        plt.plot(arange(i1*n3,i1*n3+n3),h0[i1*n3:i1*n3+n3,24]/1000.)
        plt.ylim(0,6.5)
        plt.colorbar()
        if i1==2:
            plt.xlabel("Scan")
        plt.ylabel("Height(km)")
    plt.savefig('orb'+f.split('.')[-3]+'.png')
    plt.show()
    dnCoeff=array([-0.01257341, -0.00933038])
    
    for i1,j1 in zip(a[0][b],a[1][b]):
        if abs(j1-24)<5 and bzd[i1,j1]<bcf[i1,j1]-4:
            zKuL.append(zKu[i1,j1,bzd[i1,j1]-60:bzd[i1,j1]+4])
            zKaL.append(zKa[i1,j1-12,bzd[i1,j1]-60:bzd[i1,j1]+4])
            hL.append(htop[i1,j1])
            zKu1=zKu[i1,j1,:]
            zKa1=zKa[i1,j1-12,:]
            imu=3
            dr=0.125
            eps=0.95
            dnst=0.2

            dn1d,dm1d,rrate1d,\
                zkuc,zkasim,epst,piaKu,piaKa =\
                    cmb.iter_profcv(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],\
                                    bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,2,dnCoeff,dnst)
            

            r1L.append([rrate1d[bcf[i1,j1]],dprsfcRate[i1,j1],\
                        pType[i1,j1],dm1d[bcf[i1,j1]],dn1d[bcf[i1,j1]],\
                        piaKu,piaSRT[i1,j1],epsDPR[i1,j1,bcf[i1,j1]],dnst])


r1L=array(r1L)
