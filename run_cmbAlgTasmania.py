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
fs=glob.glob("SP/2A*.HDF5")
icount=0
iplot=1
zKuL=[]
zKaL=[]
hL=[]
r1L=[]
ibbL=[]
rd1L=[]
for f in sorted(fs):
    fh=Dataset(f)
    lat=fh['NS/Latitude'][:,24]
    lon=fh['NS/Longitude'][:,24]
    a=nonzero(lat<-43)
    b=nonzero(lon[a]>100)
    n1=a[0][b][0]
    n2=a[0][b][-1]
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
    a=nonzero(pType==1)
    b=nonzero((h0[a]-1200)*(h0[a]-4500)<0)
    icount+=len(a[0][b])
    if len(a[0])>0:
        print(icount,h0[a].min(),h0[a].max())
    c=nonzero(a[1][b]==24)
    for i1,j1 in zip(a[0][b],a[1][b]):
        if abs(j1-24)<5 and bzd[i1,j1]<bcf[j1,j1]-4:
            zKuL.append(zKu[i1,j1,bzd[i1,j1]-60:bzd[i1,j1]+4])
            zKaL.append(zKa[i1,j1-12,bzd[i1,j1]-60:bzd[i1,j1]+4])
            hL.append(htop[i1,j1])
            zKu1=zKu[i1,j1,:]
            zKa1=zKa[i1,j1-12,:]
            imu=3
            dr=0.125
            eps=0.95
            dnst=-0.35
            if binBB[i1,j1]<=1:
                dn1d,dm1d,rrate1d,\
                    zkuc,zkasim,epst,piaKu,piaKa =\
                        cmb.iter_profcv(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,1,dnst)
                ibbL.append(0)
            else:
                bbb=binBB[i1,j1]+2
                bbt=binBBT[i1,j1]
                bb=binBB[i1,j1]
                dnst=0.0
                dnCoeff=array([ 0.01608098, -0.82884294])
                #dnCoeff=array([-0.00570364,  0.13319214])
                #dnst=-0.35
                dnp=random.randn(176)*0.025
                dn1d,dm1d,rrate1d,\
                    zkuc,zkasim,epst,piaKu,piaKa =\
                        cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,dnst,dnCoeff,dnp)
                dn1d1,dm1d1,rrate1d1,\
                    zkuc1,zkasim1,epst1,piaKu1,piaKa1 =\
                        cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,dnst-0.1,dnCoeff,dnp)
                
                ddpia=-((piaKa1-piaKu1)-(piaKa-piaKu))/0.1
                dpiaSRT=5*piaSRT[i1,j1]
                
                if reliabFlag[i1,j1]==1 or reliabFlag[i1,j1]==2:
                    dnst=(dpiaSRT-(piaKa-piaKu))*ddpia/(ddpia**2+1)
                    dn1d,dm1d,rrate1d,\
                        zkuc,zkasim,epst,piaKu,piaKa =\
                        cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,dnst,dnCoeff,dnp)
                    rd1L.append([rrate1d[bcf[i1,j1]],dprsfcRate[i1,j1],\
                                 pType[i1,j1],dm1d[bcf[i1,j1]],dn1d[bcf[i1,j1]],\
                                 piaKu,piaSRT[i1,j1],epsDPR[i1,j1,bcf[i1,j1]],dnst])
                    
               
                

                dnCoeff=array([-0.01257341, -0.00933038])
                dnst=0.0
                nEns=50
                #piaKuEns=[]
                #for k in range(nEns):
                #    dnp=random.randn(176)*0.025#+random.rand()*0.25
                dn1d,dm1d,rrate1d,\
                    zkuc,zkasim,epst,piaKu,piaKa =\
                        cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,dnst,dnCoeff,dnp)
                dn1d1,dm1d1,rrate1d1,\
                    zkuc1,zkasim1,epst1,piaKu1,piaKa1 =\
                        cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,dnst-0.1,dnCoeff,dnp)

                ddpia=-((piaKa1-piaKu1)-(piaKa-piaKu))/0.1
                dpiaSRT=5*piaSRT[i1,j1]
                ibbL.append(1)
                if reliabFlag[i1,j1]==1 or reliabFlag[i1,j1]==2:
                    piaKu0=piaKu
                    dnst=0.1*(dpiaSRT-(piaKa-piaKu))*ddpia/(ddpia**2+1)
                    dn1d,dm1d,rrate1d,\
                    zkuc,zkasim,epst,piaKu,piaKa =\
                        cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                        bsfc[i1,j1],zKu1,zKa1,dr,eps,imu,dnst,dnCoeff,dnp)
                    r1L.append([rrate1d[bcf[i1,j1]],dprsfcRate[i1,j1],\
                                pType[i1,j1],dm1d[bcf[i1,j1]],dn1d[bcf[i1,j1]],\
                                piaKa-piaKu,dpiaSRT,piaKu0,epsDPR[i1,j1,bcf[i1,j1]],dnst])
#                    if piaKu>0.5:
#                        stop
                
    if len(c[0])>4000:
        iplot=1
        print(c[0])
    else:
        iplot=0
        
    if iplot==1:
        nx=24
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
        plt.show()

r1L=array(r1L)
rd1L=array(rd1L)
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.loglog(r1L[:,0],rd1L[:,0],'*')
ax.plot([0.1,30],[0.1,30])
plt.xlabel('SurfaceRain=f(Protat R-Dm)')
plt.ylabel('SurfaceRain=f(DPR stratiform R-Dm)')
plt.title('Retrieved Surface Rain Rate')
plt.savefig('surfaceRain.png')
plt.figure()
plt.scatter(rd1L[:,0],rd1L[:,-1])
plt.xlabel("SurfaceRain=f(DPR stratiform R-Dm)")
plt.ylabel("DPR epsilon")
plt.show()

stop
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


    
for (i1,j1) in zip(a[0],a[1]):
    try:
        if j1>=12 and j1<37:
            zKu1=zKu[i1,j1,:]
            zKa1=zKa[i1,j1-12,:]
            if bzd[i1,j1]<bcf[i1,j1]:# and (reliabFlag[i1,j1]==1 or reliabFlag[i1,j1]==2):
                it+=1
                imu=3
                dr=0.125
                eps=0.5
                if pType[i1,j1]==2 or pType[i1,j1]==3:
                    dn1d,dm1d,rrate1d,\
                        zkuc,zkasim,epst,piaKu,piaKa =\
                            cmb.iter_profcv(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],\
                                            bsfc[i1,j1],zKu1,zKa1,dr,eps,imu)
                if pType[i1,j1]==1:
                    if binBB[i1,j1]<=1:
                        dn1d,dm1d,rrate1d,\
                            zkuc,zkasim,epst,piaKu,piaKa =\
                                cmb.iter_profcv(bst[i1,j1],bzd[i1,j1],bcf[i1,j1],\
                                                bsfc[i1,j1],zKu1,zKa1,dr,eps,imu)
                    else:
                        bbb=binBB[i1,j1]+2
                        bbt=binBBT[i1,j1]
                        bb=binBB[i1,j1]
                        dn1d,dm1d,rrate1d,\
                            zkuc,zkasim,epst,piaKu,piaKa =\
                                cmb.iter_profst(bst[i1,j1],bzd[i1,j1],bb,bbt,bbb,bcf[i1,j1],\
                                                bsfc[i1,j1],zKu1,zKa1,dr,eps,imu)
                        
                r1L.append([rrate1d[bcf[i1,j1]],dprsfcRate[i1,j1],\
                            pType[i1,j1],dm1d[bcf[i1,j1]],dn1d[bcf[i1,j1]],\
                            piaKu,piaSRT[i1,j1]])

    except:
        print('exception',i1,j1)

rrate1dL=[]
zKaSimL=[]
print('here')
for i1 in range(138):
    rrate1d=zeros(176)
    j1=24
    zKu1=zKu[i1,j1,:]
    zKa1=zKa[i1,j1-12,:]
    zkasim=zeros(176)-99.9
    eps=0.95
    #print(i1)
    if bzd[i1,24]>bcf[i1,24]:
        bzd[i1,24]=bcf[i1,24]+1
    try:
        if pType[i1,24]==2:
            a=1
            if bzd[i1,24]>bcf[i1,24]:
                bzd[i1,24]=bcf[i1,24]+1
            #print(bst[i1,24],bzd[i1,24],bcf[i1,24])
            dn1d,dm1d,rrate1d,\
                zkuc,zkasim,epst,piaKu,piaKa =\
                    cmb.iter_profcv(bst[i1,24],bzd[i1,24],bcf[i1,24],\
                                    bsfc[i1,24],zKu1,zKa1,dr,eps,imu)
            #print(i1)
        if pType[i1,24]==1:
            if binBB[i1,j1]<=1:
                a=1
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
                        cmb.iter_profst(bst[i1,24],bzd[i1,24],bb,bbt,\
                                        bbb,bcf[i1,24],\
                                        bsfc[i1,24],zKu1,zKa1,dr,eps,imu)
    except:
        print(i1)
        pass
        
    rrate1dL.append(rrate1d)
    zKaSimL.append(zkasim)


rrate2d=array(rrate1dL)
rrate2d=ma.array(rrate2d,mask=rrate2d<0.001)
plt.figure()
plt.pcolormesh(rrate2d[:,::-1].T,cmap='jet',vmax=30)
plt.ylim(0,60)
#plt.xlim(200,400)
plt.colorbar()


zKaSimL=array(zKaSimL)
zKaSimL=ma.array(zKaSimL,mask=zKaSimL<0.001)
plt.figure()
plt.subplot(211)
plt.pcolormesh(zKaSimL[:,::-1].T,cmap='jet',vmin=5,vmax=40)
plt.ylim(0,60)
#plt.xlim(0,700)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(zKa[:,12,::-1].T,cmap='jet',vmin=5,vmax=40)
plt.ylim(0,60)
#plt.xlim(0,700)
plt.colorbar()
plt.show()

