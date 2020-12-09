fname='2A.GPM.DPR.V8-20180723.20180118-S125407-E142641.022105.V06A.HDF5'
fname='2A.GPM.DPR.V8-20180723.20180215-S134018-E151250.022541.V06A.HDF5'

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
fst=[]
#for i in range(30):
#    fs=glob.glob("/gpmdata/2018/08/%2.2i/radar/2A.GPM.DPR.V8*HDF5"%(i+1))
#    fst.extend(sorted(fs))
fst=glob.glob("../monthly/SEAsia/2A*DPR*HDF5")
fs=glob.glob("../monthly/SEAsiaCS/2A*DPR*HDF5")
fst=sorted(fst)
fst.extend(sorted(fs))
fstKu=glob.glob("../ORO/DPR-CS/2A-CS-CON*DPR*HDF5")
fstKu=glob.glob("MediCane/2A*Ku*HDF5")
fstKu=glob.glob("../monthly/SEAsiaCS/2A*DPR*HDF5")
fstKu=glob.glob("../monthly/SEAsia/2A*DPR*HDF5")
fstKu=glob.glob("../cmbv6x/SP/2A*DPR*HDF5")
print(len(fst))
#stop
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
fst=fst[:]

deltaSRT_2d=zeros((360,130,2,3),float)
sfcRainRate_2d=zeros((360,130,2,3),float)
count=zeros((360,130,3),float)

dataL=[]
for (ifile,f) in zip(range(0,len(fstKu)),sorted(fstKu)[0:]):
    f=f.replace("Ku","DPR")
    fh=Dataset(f)
    lat=fh['NS/Latitude'][:,24]
    lon=fh['NS/Longitude'][:,24]
    n1=0
    n2=lon.shape[0]
    rrate3d=zeros((n2-n1+1,24,176),float)
    rratesfc=zeros((n2-n1+1,24),float)
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
    sfcType=fh['NS/PRE/landSurfaceType'][:,:]
    fh.close()
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
    iplot=0
    if iplot==1:
        plt.figure(figsize=(12,8))
        for i1 in range(1,3):
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
                plt.savefig('orbKu'+f.split('.')[-3]+'.png')
        plt.figure(figsize=(12,8))
        for i1 in range(1,3):
            plt.suptitle('Orbit '+f.split('.')[-3]+' %i'%ifile)
            plt.subplot(3,1,i1+1)
            plt.pcolormesh(arange(i1*n3,i1*n3+n3),arange(176)*0.125,\
                           zKa[i1*n3:i1*n3+n3,12,::-1].T,vmin=0,vmax=40,cmap="jet")
            plt.plot(arange(i1*n3,i1*n3+n3),h0[i1*n3:i1*n3+n3,24]/1000.)
            plt.ylim(0,12.5)
            plt.colorbar()
            if i1==1:
                plt.xlabel("Scan")
                plt.ylabel("Height(km)")
                plt.savefig('orbKa'+f.split('.')[-3]+'.png')
        plt.show()
    ifile+=1
    
    dnCoeff=array([-0.01257341, -0.00933038])
    j1=24
    zKaSim=zeros((3*n3+3,176),float)
    zKuSim=zeros((3*n3+3,176),float)
    rrate2d=zeros((3*n3+3,176),float)
    a=nonzero(pType>0)
    b=nonzero(abs(a[1]-24)<12)
    for i1,j1 in zip(a[0][b],a[1][b]):
        zKu1=zKu[i1,j1,:]
        zKa1=zKa[i1,j1-12,:]
        imu=3
        dr=0.125
        eps=0.95
        dnst=0.2
        if pType[i1,j1]==0:
            continue
        if(reliabFlag[i1,j1]==2):
            reliabFlag[i1,j1]=1
        if pType[i1,j1]==2:
            reliabFlag[i1,j1]=0
        dsrtPIA=piaSRT[i1,j1]*5
        dn=0
        eps=1.0
        if bzd[i1,j1]<bcf[i1,j1]+40:
            if pType[i1,j1]==2 and abs(lat[i1,j1])<45:
                if h0[i1,j1]>2000:
                    dn=0.2*(h0[i1,j1]-2000.0)/2000.0
            dn1d,dm1d,rrate1d,zkuc,zkasim,\
                epst,piaku,\
                piaka,dt1,dt2 = cmb.prof1d(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],bsfc[i1,j1],\
                                           binBB[i1,j1],binBBT[i1,j1],zKu1,zKa1,pType[i1,j1],\
                                           dr,eps,imu,dnCoeff,dn,\
                                           dsrtPIA,reliabFlag[i1,j1],piaSRT[i1,j1],\
                                           -1)#,[n1d])
            
            
            dataL.append([rrate1d,dprsfcRate[i1,j1],bst[i1,j1],bcf[i1,j1],lon[i1,j1],lat[i1,j1],\
                          reliabFlag[i1,j1],dsrtPIA,piaka,piaku,zKu1,pType[i1,j1],sfcType[i1,j1]])
            #if pType[i1,j1]==1 and reliabFlag[i1,j1]!=1:
            #    print('here')
            ix=int(lon[i1,j1]+180)
            jx=int(lat[i1,j1]+65)
            rrate3d[i1,j1-12,:]=rrate1d
            rratesfc[i1,j1-12]=rrate1d[bcf[i1,j1]]
            if ix>=0 and ix<360 and jx>=0 and jx<130:
                if dprsfcRate[i1,j1]>-0.1:
                    #print(reliabFlag[i1,j1],pType[i1,j1],dprsfcRate[i1,j1],ix,jx)
                    if reliabFlag[i1,j1]!=1:
                        if pType[i1,j1]==1:
                            sfcRainRate_2d[ix,jx,0,0]+=dprsfcRate[i1,j1]
                            sfcRainRate_2d[ix,jx,1,0]+=rrate1d[bcf[i1,j1]]
                            count[ix,jx,0]+=1
                            deltaSRT_2d[ix,jx,0,0]+=dsrtPIA
                            deltaSRT_2d[ix,jx,1,0]+=(piaka-piaku)#dsrtPIA
                        if pType[i1,j1]==2:
                            sfcRainRate_2d[ix,jx,0,1]+=dprsfcRate[i1,j1]
                            sfcRainRate_2d[ix,jx,1,1]+=rrate1d[bcf[i1,j1]]
                            count[ix,jx,1]+=1
                            deltaSRT_2d[ix,jx,0,1]+=dsrtPIA
                            deltaSRT_2d[ix,jx,1,1]+=(piaka-piaku)#dsrtPIA
                    if reliabFlag[i1,j1]==1:
                        if pType[i1,j1]==1:
                            sfcRainRate_2d[ix,jx,0,2]+=dprsfcRate[i1,j1]
                            sfcRainRate_2d[ix,jx,1,2]+=rrate1d[bcf[i1,j1]]
                            count[ix,jx,2]+=1
                            deltaSRT_2d[ix,jx,0,2]+=dsrtPIA
                            deltaSRT_2d[ix,jx,1,2]+=(piaka-piaku)#dsrtPIA
                    
    #stop
import pickle

pickle.dump(dataL,open('sp_noG_Data.pklz','wb'))
