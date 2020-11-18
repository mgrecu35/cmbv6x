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
#fs=glob.glob("SP/2A*SP*HDF5")
icount=0
iplot=1
zKuL=[]
zKaL=[]
hL=[]
r1L=[]
ibbL=[]
rd1L=[]
#for f in sorted(fs)[0:1]:
iplot=0
ifile=0
dFL=[]
DPR_RDmCoeff=[ 0.01608098, -0.82884294]
for (ifile,f) in zip(range(0,53),sorted(fs)[0:53]):
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
        #plt.suptitle(f)
    print(ifile,f)
    if iplot==1:
        plt.figure(figsize=(12,8))
        for i1 in range(3):
            plt.suptitle('Orbit '+f.split('.')[-3]+' %i'%ifile)
            plt.subplot(3,1,i1+1)
            plt.pcolormesh(arange(i1*n3,i1*n3+n3),arange(176)*0.125,\
                           zKu[i1*n3:i1*n3+n3,24,::-1].T,vmin=0,vmax=45,cmap="jet")
            plt.plot(arange(i1*n3,i1*n3+n3),h0[i1*n3:i1*n3+n3,24]/1000.)
            plt.ylim(0,12.5)
            plt.colorbar()
            if i1==1:
                plt.xlabel("Scan")
                plt.ylabel("Height(km)")
                plt.savefig('orb'+f.split('.')[-3]+'.png')
        plt.show()
    ifile+=1
    dnCoeff=array([-0.01257341, -0.00933038])
    j1=24
    zKaSim=zeros((3*n3,176),float)
    zKuSim=zeros((3*n3,176),float)
    rrate2d=zeros((3*n3,176),float)
    for i1 in range(3*n3):
        zKu1=zKu[i1,j1,:]
        zKa1=zKa[i1,j1-12,:]
        imu=3
        dr=0.125
        eps=0.95
        dnst=0.2
        if pType[i1,j1]==0:
            continue
        if pType[i1,j1]==2:
            dnCoeff=array([ 0.01608098, -0.82884294])
            dnCoeff=array([ 0.00827692, -0.19116723])
            dn=0.2
        else:
            dnCoeff=array([-0.00570364,  0.13319214])
            dnCoeff=array([ 0.01608098, -0.82884294])
            dnCoeff=array([-0.01238408, -0.01593891])
            dn=0.0+0.2
        if(reliabFlag[i1,j1]==2):
            reliabFlag[i1,j1]=1
        if pType[i1,j1]==2:
            reliabFlag[i1,j1]=0
        dsrtPIA=piaSRT[i1,j1]*5
        if bzd[i1,j1]<bcf[i1,j1]+40:
            dn1d,dm1d,rrate1d,zkuc,zkasim,\
                epst,piaku,piaka = cmb.prof1d(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],bsfc[i1,j1],\
                                              binBB[i1,j1],binBBT[i1,j1],zKu1,zKa1,pType[i1,j1],dr,eps,imu,dnCoeff,dn,\
                                              dsrtPIA,reliabFlag[i1,j1])#,[n1d])
            
            a=nonzero(zKa1[bzd[i1,j1]+4:bcf[i1,j1]]>10)
            b=nonzero(zKu1[bzd[i1,j1]+4:bcf[i1,j1]][a]>10)
            if len(b[0])>0 and pType[i1,j1]==1:
                dF=(((zKa1[bzd[i1,j1]+4:bcf[i1,j1]][a][b]-zkasim[bzd[i1,j1]+4:bcf[i1,j1]][a][b])**2).mean())**0.5
                #print(dF)
                dFL.append(dF)
            zKaSim[i1,:]=zkasim
            zKuSim[i1,:]=zkuc
            rrate2d[i1,:]=rrate1d
        #if reliabFlag[i1,j1]==1:
            r1L.append([rrate1d[bcf[i1,j1]],dprsfcRate[i1,j1],dm1d[bcf[i1,j1]],dn1d[bcf[i1,j1]],dsrtPIA,piaka-piaku,pType[i1,j1]])

zKaSim=ma.array(zKaSim,mask=zKaSim<10)
zKuSim=ma.array(zKuSim,mask=zKuSim<10)
rrate2d=ma.array(rrate2d,mask=rrate2d<0.01)
import matplotlib
for i1 in range(0,3):
    plt.figure()
    plt.subplot(3,1,1)
    plt.pcolormesh(arange(i1*n3,i1*n3+n3),arange(176)*0.125,\
                   zKu[i1*n3:i1*n3+n3,24,::-1].T,vmin=0,vmax=50,cmap="jet")
    plt.plot(arange(i1*n3,i1*n3+n3),h0[i1*n3:i1*n3+n3,24]/1000.)
    plt.ylim(0,12.5)
    plt.colorbar()
    plt.ylabel("Height(km)")
    plt.subplot(3,1,2)
    plt.pcolormesh(arange(i1*n3,i1*n3+n3),arange(176)*0.125,\
                   zKuSim[i1*n3:i1*n3+n3,::-1].T,vmin=0,vmax=50,cmap="jet")
    plt.plot(arange(i1*n3,i1*n3+n3),h0[i1*n3:i1*n3+n3,24]/1000.)
    plt.ylim(0,12.5)
    plt.colorbar()
    plt.subplot(3,1,3)
    plt.pcolormesh(arange(i1*n3,i1*n3+n3),arange(176)*0.125,\
                   rrate2d[i1*n3:i1*n3+n3,::-1].T,vmin=0.1,vmax=100,norm=matplotlib.colors.LogNorm(),cmap="jet")
    plt.ylim(0,12.5)
    plt.colorbar()
