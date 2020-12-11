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
import pickle
import time

iread=1

import glob
fs=sorted(glob.glob("pklzs/2018AugSt*5000_*pklz"))
tfort=0
tpyth=0
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
    
    print('pklzs in done',len(dataL),time.time()-t1)
    iproft=0
    iprofBB=0
    for d1 in dataL:
        zKu,zKa,bst,bcf,bzd,binBB,binBBT,lon,lat,\
            reliabFlag,piaSRT,pType,sfcType,dprsfcRate,bsfc,\
            relFlagKu,piaSRTKu,cmbsfcRain=d1
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
        if bzd+6<bcf:
            #binBBT=0
            #binBB=bzd-1
            #binBBT=binBB-3
            dn=-0.1
            t1st=time.time()
            dn1d,dm1d,rrate1d,zkuc,zkasim,\
                epst,piaku,\
                piaka,dt1,dt2 = cmb.prof1d(bst,bzd,bcf,bsfc,\
                                           binBB,binBBT,zKu,zKa,pType,\
                                           dr,eps,imu,dnCoeff,dn,
                                           dsrtPIA,reliabFlag,-99.9,\
                                           -1)
            dn1d,dm1d,rrate1d_noad,zkuc,zkasim_noad,\
                epst,piaku_noad,\
                piaka_noad,dt1,dt2 = cmb.prof1d_noad(bst,bzd,bcf,bsfc,\
                                                     binBB,binBBT,zKu,zKa,pType,\
                                                     dr,eps,imu,dnCoeff,dn,
                                                     dsrtPIA,reliabFlag,-99.9,\
                                                     -1)
            iproft+=1
            t1st=time.time()-t1st
            tfort+=t1st
            if binBB>0 and binBB<-170 and rrate1d[bcf]>0.1 and zKu[binBB]>10:
                iprofBB+=1
                bbb=binBB+2
                dnp=zeros((176),float)
                dn1dd,dm1dd,rrate1dd,zkucd,zkasimd,\
                    epstd,piakud,piakad,\
                    dzdnd,dpiadnd = cmb.iter_profst(bst,bzd,binBB,binBBT,bbb,bcf,bsfc,\
                                                    zKu,zKa,dr,eps,imu,dn,dnCoeff,dnp)
                dzL=[]
                t1st=time.time()
                nz=bcf+1-bbb
                jacob=zeros((nz+1,nz))
                dy=zeros((nz+1),float)
                for k in range(bbb,bcf+1):
                    dnp[bbb:bcf+1]=0
                    dnp[k]=-0.1
                    dn1dp,dm1dp,rrate1dp,zkucp,zkasimp,\
                    epstp,piakup,piakap,\
                    dzdnp,dpiadnp = cmb.iter_profst(bst,bzd,binBB,binBBT,bbb,bcf,bsfc,\
                                                    zKu,zKa,dr,eps,imu,dn,dnCoeff,dnp)
                    dzL.append((zkasimp[bcf]-zkasimd[bcf])/(-0.1))
                    jacob[:-1,k-bbb]=(zkasimp[bbb:bcf+1]-zkasimd[bbb:bcf+1])/(-0.1)
                    jacob[-1,k-bbb]=((piakap-piakup)-(piakad-piakud))/(-0.1)
                    if zKu[k]>10 and zKa[k]>10:
                        dy[k-bbb]=zKa[k]-zkasimd[k]
                if reliabFlag==1:
                    dy[-1]=dsrtPIA-(piakad-piakud)

                tmp=dot(jacob.T,jacob)
                graddy=dot(jacob.T,dy)
                tmp=tmp+eye(nz)*4
                dnsol=linalg.solve(tmp,graddy)
                dnsol[dnsol<-1.5]=-1.5
                dnsol[dnsol>1.5]=1.5
                dnp[bbb:bcf+1]=dnsol
                dn1d,dm1d,rrate1d,zkuc,zkasim,\
                    epst,piaku,piaka,\
                    dzdn,dpiadn = cmb.iter_profst(bst,bzd,binBB,binBBT,bbb,bcf,bsfc,\
                                                    zKu,zKa,dr,eps,imu,dn,dnCoeff,dnp)
                #print(dnsol)
                t1st=time.time()-t1st
                tpyth+=t1st
            if reliabFlag==1:
                r1L.append([rrate1d[bcf],dprsfcRate,piaka-piaku,dsrtPIA,sfcType,binBB,lat,cmbsfcRain,\
                            rrate1d_noad[bcf],piaka_noad-piaku_noad,zKa,zkasim,zkasim_noad,bzd,binBB,bcf])
                bzdL.append(bzd)
                zku1L.append(zKu[bzd-56:bzd+6])
                zka1L.append(zKa[bzd-56:bzd+6])
            else:
                r1Lun.append([rrate1d[bcf],dprsfcRate,piaka-piaku,dsrtPIA,sfcType,binBB,lat,cmbsfcRain,
                              rrate1d_noad[bcf],piaka_noad-piaku_noad,zKa,zkasim,zkasim_noad,bzd,binBB,bcf])


    pickle.dump([r1L,r1Lun],open('retrs/stratRetrievalsAug2018St_%3.3i.pklz'%ifile,'wb'))
    print("tfort=%f tpyth=%f #profs=(%i,%i)"%(tfort,tpyth,iproft,iprofBB))
    ifile+=1
