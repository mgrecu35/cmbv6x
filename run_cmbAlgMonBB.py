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
            dt1=0
            dt2=0
            dnp=zeros((176),float)
            if binBB>0:
                bbb=binBB+2
                dn1d,dm1d,rrate1d,zkuc,zkasim,\
                    epst,piaku,piaka,\
                    dzdn,dpiadn = cmb.iter_profst(bst,bzd,binBB,binBBT,bbb,bcf,bsfc,\
                                                  zKu,zKa,dr,eps,imu,dn,dnCoeff,dnp)
                dn1d1,dm1d1,rrate1d1,zkuc,zkasim1,\
                    epst1,piaku1,piaka1,\
                    dzdn1,dpiadn1 = cmb.iter_profst(bst,bzd,binBB,binBBT,bbb,bcf,bsfc,\
                                                  zKu,zKa,dr,eps,imu,dn+0.1,dnCoeff,dnp)
                gradZS=1
                gradZe=0.
                for k in range(bbb+3,bcf+1):
                    if zKu[k]>10.0 and zKa[k]>10.0:
                        gradZS=gradZS+0.125*((zkasim1[k]-zkasim[k])/0.1)**2
                        gradZe=gradZe+0.25*(zkasim1[k]-zkasim[k])/0.1*0.25*(zKa[k]-zkasim[k])
                        
                        
                if reliabFlag==1:
                    dpiadn=(piaka1-piaku1-piaka+piaku)/(0.1)
                    gradZS=gradZS+dpiadn**2
                    gradZe=gradZe-dpiadn*(piaka-piaku-dsrtPIA)
                    
                ddn=0.95*gradZe/gradZS
                dn1d,dm1d,rrate1d,zkuc,zkasim,\
                    epst,piakuf,piakaf,\
                    dzdn,dpiadn = cmb.iter_profst(bst,bzd,binBB,binBBT,bbb,bcf,bsfc,\
                                                  zKu,zKa,dr,eps,imu,dn+ddn,dnCoeff,dnp)
                #print(piaka-piaku,piaka1-piaku1,dsrtPIA,piakaf-piakuf,eps)
            else:
                dn1d,dm1d,rrate1d,zkuc,zkasim,\
                    epst,piaku,piaka,\
                    dzdn,dpiadn,dt1,dt2 = cmb.iter_profst_nobb(bst,bzd,bcf,bsfc,\
                                                               zKu,zKa,dr,eps,imu,pType,dnCoeff,dn,dnp)
                
            if reliabFlag==1:
                r1L.append([rrate1d[bcf],dprsfcRate,piaka-piaku,dsrtPIA,sfcType,binBB,lat,cmbsfcRain,\
                            rrate1d[bcf],piaka-piaku,zKa,zkasim,zkasim,bzd,binBB,bcf])
                bzdL.append(bzd)
                zku1L.append(zKu[bzd-56:bzd+6])
                zka1L.append(zKa[bzd-56:bzd+6])
            else:
                r1Lun.append([rrate1d[bcf],dprsfcRate,piaka-piaku,dsrtPIA,sfcType,binBB,lat,cmbsfcRain,
                              rrate1d[bcf],piaka-piaku,zKa,zkasim,zkasim,bzd,binBB,bcf])


    pickle.dump([r1L,r1Lun],open('retrs/stratRetrievalsAug2018St_%3.3i.pklz'%ifile,'wb'))
    print("tfort=%f tpyth=%f #profs=(%i,%i)"%(tfort,tpyth,iproft,iprofBB))
    ifile+=1
