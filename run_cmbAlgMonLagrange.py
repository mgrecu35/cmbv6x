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

import xarray as xr
zkuxr=xr.DataArray(cmb.tablep2.zkusj)
zkaxr=xr.DataArray(cmb.tablep2.zkasj)
rjxr=xr.DataArray(cmb.tablep2.rj)
dmjr=xr.DataArray(cmb.tablep2.dmj)
attkujr=xr.DataArray(cmb.tablep2.attkuj)
attkajr=xr.DataArray(cmb.tablep2.attkaj)

zku=cmb.tablep2.zkusj
rj=cmb.tablep2.rj
dm=cmb.tablep2.dmj
attku=cmb.tablep2.attkuj
attka=cmb.tablep2.attkaj
import numpy as np
attKaCoeff=np.polyfit(zkuxr[110:169],log(attkajr[110:169]),2)
attKuCoeff=np.polyfit(zkuxr[110:189],log(attkujr[110:189]),2)
pred=exp(np.polyval(attKaCoeff,zkuxr[:110]))
predKu=exp(np.polyval(attKuCoeff,zkuxr[:120]))
plt.semilogx(predKu[10:120],zkuxr[10:120])
plt.semilogx(attkujr[10:189],zkuxr[10:189])

attkajr[:110]=pred
attkujr[:120]=predKu
tablexr=xr.Dataset({"zKu":zkuxr,"Dm":dmjr,"attKu":attkujr,"attKa":attkajr,"rainRate":rjxr})
table_cmb={"zKu":zkuxr,"zKa":zkaxr,"Dm":dmjr,"attKu":attkujr,"attKa":attkajr,"rainRate":rjxr}
import pandas as pd
#df = pd.DataFrame(data=table_cmb)
#df.to_csv("cmbTables.csv")

#from rbf.interpolate import RBFInterpolant
#zkuint = RBFInterpolant(log(rjxr), zkuxr, sigma=0.1, phi='ga', order=1)
# create the interpolation points, and evaluate the interpolant
#zku_itp = zkuint(log(rjxr)) 
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
import pickle
import time

iread=1

import glob
fs=sorted(glob.glob("pklzs/2018AugSt*5000_*pklz"))
tfort=0
tpyth=0
from procOrb import *
imax=0
ic=0
for f1 in []:#fs:
    r1L=[]
    r1Lun=[]
    bzdL=[]
    zku1L=[]
    zka1L=[]
    iread=1
    print('reading pkzls %s'%f1)
    t1=time.time()
    dataL=pickle.load(open(f1,'rb'))
    
    print('pklzs in done',len(dataL),time.time()-t1)
    iproft=0
    iprofBB=0
    contigPix,nc,zku3d,zka3d=procorb(dataL)
    max1=0
    for i in range(nc):
        a1=nonzero(contigPix==i+1)
        if len(a1[0])>100:
            print(a1[0][0],a1[0][-1])
            ind=argmax(contigPix[a1[0][0]:a1[0][-1],:].sum(axis=0))
            if a1[0][-1]-a1[0][0]>50:
                plt.figure(figsize=(12,8))
                plt.subplot(211)
                plt.pcolormesh(zku3d[a1[0][0]:a1[0][-1],ind,::-1].T,\
                               cmap='jet',vmin=0,vmax=50)
                plt.ylim(0,80)
                plt.subplot(212)
                plt.pcolormesh(zka3d[a1[0][0]:a1[0][-1],ind,::-1].T,\
                               cmap='jet',vmin=0,vmax=50)
                plt.ylim(0,80)
                ic+=1
            if a1[0][-1]-a1[0][0]>max1:
                max1=a1[0][-1]-a1[0][0]
    if max1>100:
        imax+=1
    if imax==30 or ic>100:
        stop
for f1 in fs:
    r1L=[]
    r1Lun=[]
    bzdL=[]
    zku1L=[]
    zka1L=[]
    iread=1
    print('reading pkzls %s'%f1)
    t1=time.time()
    dataL=pickle.load(open(f1,'rb'))
    for d1 in dataL:
        zKu,zKa,bst,bcf,bzd,binBB,binBBT,lon,lat,\
            reliabFlag,piaSRT,pType,sfcType,dprsfcRate,bsfc,\
            relFlagKu,piaSRTKu,cmbsfcRain,piaHY,i1,j1,pType=d1
        #print(reliabFlag)
        dnCoeff=array([-0.01257341, -0.00933038])
        imu=3
        dr=0.125
        if(reliabFlag==2):
            reliabFlag=1
        if pType==2:
            reliabFlag=0
        dsrtPIA=piaSRT*5
        dn=0
        eps=1.0
        if reliabFlag==1:
            reliabFlag=1
        icount=0
        if bzd+6<bcf and pType==1:
            #binBBT=0
            #binBB=bzd-1
            #binBBT=binBB-3
            dn=-0.1
            t1st=time.time()
            dt1=0
            dt2=0
            dnp=zeros((176),float)
            if binBB>0:
                bbb=binBB+2
                dn1d,dm1d,rrate1d,zkuc,zkasim,\
                    epst,piaku,piaka,\
                    dzdn,dpiadn,\
                    piaKuS,piaKaS = cmb.iter_profst(bst,bzd,binBB,binBBT,\
                                                bbb,bcf,bsfc,\
                                                zKu,zKa,dr,eps,imu,dn,\
                                                dnCoeff,dnp)
                nc=bsfc-bcf
                wzku=0.2
                wzka=0.1
                wpia=1
                nens=60
                dpiaRet=piaka-piaku
                if bcf>bbb+1:
                    rrate_out,dn_out,zkusimE,zkasimE,\
                        zku_out,zka_out,rrEns,yEns,xEns,dy,pia_out = \
                        cmb.rainprofstg(zKu[bbb:bcf+1],zKa[bbb:bcf+1],\
                                        dsrtPIA,piaKuS,piaKaS,\
                                        reliabFlag,nc,dr,wzku,wzka,wpia,\
                                        rrate1d[bbb:bcf+1],dn1d[bbb:bcf+1],nens)
                    xEns=xEns.T
                    yEns=yEns.T
                    nx=xEns.shape[0]
                    covXY=cov(xEns,yEns)[:nx,nx:]
                    covYY=cov(yEns)
                    ny=yEns.shape[0]
                    R=4.
                    rrate1dF=rrate1d.copy()
                    kGain=dot(covXY,linalg.inv(covYY+R*eye(ny)))
                    xRet=xEns.mean(axis=1)+dot(kGain,dy)
                    #print(dot(kGain,dy))
                    #stop
                    xMin=xEns.min(axis=1)
                    xMax=xEns.max(axis=1)
                    a=nonzero(xRet<xMin)
                    xRet[a]=xMin[a]
                    a=nonzero(xRet>xMax)
                    xRet[a]=xMax[a]
                    n1=int((nx-1)/2)
                    rrate1dF[bbb:bcf+1]=rrate_out#xRet[0:n1]
                    zkasimF=zkasim.copy()
                    zkasimF[bbb:bcf+1]=zka_out
                    #stop
                    dpiaRet=xRet[-1]
                    dpiaRet=pia_out[1]-pia_out[0]
            else:
                dn1d,dm1d,rrate1d,zkuc,zkasim,\
                    epst,piaku,piaka,\
                    dzdn,dpiadn,\
                    dt1,dt2,\
                    piaKuS,piaKaS = cmb.iter_profst_nobb(bst,bzd,bcf,bsfc,\
                                                         zKu,zKa,dr,eps,imu,\
                                                         pType,dnCoeff,dn,dnp)
                
                bbb=bzd+1
                nc=bsfc-bcf
                wzku=0.2
                wzka=0.1
                wpia=1
                nens=60
                dpiaRet=piaka-piaku
                if bcf>bbb+1:
                    rrate_out,dn_out,zkusimE,zkasimE,\
                        zku_out,zka_out,rrEns,yEns,xEns,dy,pia_out = \
                        cmb.rainprofstg(zKu[bbb:bcf+1],zKa[bbb:bcf+1],\
                                        dsrtPIA,piaKuS,piaKaS,\
                                        reliabFlag,nc,dr,wzku,wzka,wpia,\
                                        rrate1d[bbb:bcf+1],dn1d[bbb:bcf+1],nens)
                    xEns=xEns.T
                    yEns=yEns.T
                    nx=xEns.shape[0]
                    covXY=cov(xEns,yEns)[:nx,nx:]
                    covYY=cov(yEns)
                    ny=yEns.shape[0]
                    R=4.
                    kGain=dot(covXY,linalg.inv(covYY+R*eye(ny)))
                    xRet=xEns.mean(axis=1)+dot(kGain,dy)
                    xMin=xEns.min(axis=1)
                    xMax=xEns.max(axis=1)
                    a=nonzero(xRet<xMin)
                    xRet[a]=xMin[a]
                    a=nonzero(xRet>xMax)
                    xRet[a]=xMax[a]
                    n1=int((nx-1)/2)
                    rrate1dF=rrate1d.copy()
                    zkasimF=zkasim.copy()
                    zkasimF[bbb:bcf+1]=zka_out
                    rrate1dF[bbb:bcf+1]=rrate_out
                    dpiaRet=xRet[-1]
                    dpiaRet=pia_out[1]-pia_out[0]
            if reliabFlag==1:
                if(zKa[bcf]>10 and zKu[bcf]>10):
                    r1L.append([rrate1d[bcf],dprsfcRate,piaka-piaku,\
                                dsrtPIA,sfcType,binBB,lat,cmbsfcRain,\
                                rrate1d[bcf],piaka-piaku,zKa,zkasim,zkasim,\
                                bzd,binBB,bcf,rrate1dF[bcf],zkasimF[bcf],\
                                dpiaRet])
                if piaka-piaku==inf:
                    stop
                bzdL.append(bzd)
                zku1L.append(zKu[bzd-56:bzd+6])
                zka1L.append(zKa[bzd-56:bzd+6])
            else:
                p=1
                if(zKa[bcf]>10 and zKu[bcf]>10):
                    r1Lun.append([rrate1d[bcf],dprsfcRate,piaka-piaku,\
                                  dsrtPIA,sfcType,binBB,lat,cmbsfcRain,
                                  rrate1d[bcf],piaka-piaku,zKa,zkasim,zkasim,\
                                  bzd,binBB,bcf,rrate1dF[bcf],zkasimF[bcf],\
                                  dpiaRet])


    pickle.dump([r1L,r1Lun],open('retrs/stratRetrievalsAug2018StE_%3.3i.pklz'%ifile,'wb'))
    #print("tfort=%f tpyth=%f #profs=(%i,%i)"%(tfort,tpyth,iproft,iprofBB))
    ifile+=1
