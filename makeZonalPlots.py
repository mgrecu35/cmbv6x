iread=1
from numpy import *
import matplotlib.pyplot as plt
import pickle
import glob

fs=sorted(glob.glob("retrs/stratRetrievalsAug2018St*"))
dx=2.5
nx=int(120/dx)
countsLand=zeros((nx,2),float)
countsOcean=zeros((nx,2),float)
rainLand=zeros((nx,2,4),float)
rainOcean=zeros((nx,2,4),float)
dpiaLand=zeros((nx,2,3),float)
dpiaOcean=zeros((nx,2,3),float)
rs=[]
rsL=[]
zsfL=[]
zsfcO=[]
for f1 in fs:
    [r1L,r1Lun]=pickle.load(open(f1,'rb'))

    for r1 in r1L:
        rate,dprrate,dpia,dpiaSRT,sfcType,binBB,lat,cmbSfcRain,\
            rate_noad,dpia_noad,zKa,zkasim,zkasim_noad,bzd,binBB,bcf=r1
        i0=int((lat+60)/dx)
        i0=max(0,i0)
        i0=min(i0,nx-1)
        if binBB>0:
            ibb=0
        else:
            ibb=1
        if sfcType==0:
            rs.append([rate,dprrate,cmbSfcRain,rate_noad])
            if zKa[bcf]>10 and binBB>0:
                zsfcO.append([zKa[bcf],zkasim[bcf],zkasim_noad[bcf]])
            countsOcean[i0,ibb]+=1
            rainOcean[i0,ibb,0]+=rate
            rainOcean[i0,ibb,3]+=rate_noad
            rainOcean[i0,ibb,1]+=dprrate
            dpiaOcean[i0,ibb,0]+=dpia
            dpiaOcean[i0,ibb,2]+=dpia_noad
            dpiaOcean[i0,ibb,1]+=dpiaSRT
            rainOcean[i0,ibb,2]+=cmbSfcRain
        else:
            rsL.append([rate,dprrate,cmbSfcRain,rate_noad])
            countsLand[i0,ibb]+=1
            rainLand[i0,ibb,0]+=rate
            rainLand[i0,ibb,3]+=rate_noad
            rainLand[i0,ibb,1]+=dprrate
            dpiaLand[i0,ibb,0]+=dpia
            dpiaLand[i0,ibb,2]+=dpia_noad
            dpiaLand[i0,ibb,1]+=dpiaSRT
            rainLand[i0,ibb,2]+=cmbSfcRain
            
a=(countsLand>0)
rmax=6
rainLand[:,:,0][a]/=countsLand[a]
rainLand[:,:,1][a]/=countsLand[a]
dpiaLand[:,:,0][a]/=countsLand[a]
dpiaLand[:,:,1][a]/=countsLand[a]
dpiaLand[:,:,2][a]/=countsLand[a]
rainLand[:,:,2][a]/=countsLand[a]
rainLand[:,:,3][a]/=countsLand[a]
a=(countsOcean>0)
rainOcean[:,:,0][a]/=countsOcean[a]
rainOcean[:,:,1][a]/=countsOcean[a]
dpiaOcean[:,:,0][a]/=countsOcean[a]
dpiaOcean[:,:,1][a]/=countsOcean[a]
dpiaOcean[:,:,2][a]/=countsOcean[a]
rainOcean[:,:,2][a]/=countsOcean[a]
rainOcean[:,:,3][a]/=countsOcean[a]
plt.suptitle('Land')
f1=plt.subplot(211)
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,0,0])
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,0,1])
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,0,2])
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,0,3],'*')
plt.ylim(0,rmax)
plt.ylabel('mm/h')
f1.axes.get_xaxis().set_visible(False)
plt.title('Bright Band')
plt.xlim(-60,60)
plt.legend(['CMB_radar','DPR','CMB_V6'])
plt.subplot(212)
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,1,0])
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,1,1])
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,1,2])
plt.plot(-60+dx/2+arange(nx)*dx,rainLand[:,1,3],'*')
plt.ylabel('mm/h')
plt.ylim(0,rmax)
plt.xlim(-60,60)
plt.title('No Bright Band')
plt.xlabel('Latitude')
plt.savefig('landStZonalRainGN_SRT.png')
plt.figure()
plt.suptitle('Ocean')
f1=plt.subplot(211)
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,0,0])
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,0,1])
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,0,2])
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,0,3],'*')
plt.ylabel('mm/h')
plt.title('Bright Band')
f1.axes.get_xaxis().set_visible(False)
plt.legend(['CMB_radar','DPR','CMB_V6'])
plt.xlim(-60,60)
plt.subplot(212)
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,1,0])
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,1,1])
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,1,2])
plt.plot(-60+dx/2+arange(nx)*dx,rainOcean[:,1,3],'*')
plt.ylabel('mm/h')
plt.title('No Bright Band')
plt.xlabel('Latitude')
plt.xlim(-60,60)
plt.savefig('oceanStZonalRainGN_SRT.png')


plt.figure()
plt.suptitle('Ocean')
f1=plt.subplot(211)
plt.plot(-60+dx/2+arange(nx)*dx,dpiaOcean[:,0,0])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaOcean[:,0,1])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaOcean[:,0,2],'*')
plt.ylabel('dB')
plt.title('Bright Band')
f1.axes.get_xaxis().set_visible(False)
plt.legend(['CMB_radar','DSRT'])
plt.xlim(-60,60)
plt.subplot(212)
plt.plot(-60+dx/2+arange(nx)*dx,dpiaOcean[:,1,0])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaOcean[:,1,1])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaOcean[:,1,2],'*')
plt.ylabel('dB')
plt.title('No Bright Band')
plt.xlabel('Latitude')
plt.xlim(-60,60)
plt.savefig('oceanStZonalDPIAGN_SRT.png')


plt.figure()
plt.suptitle('Land')
f1=plt.subplot(211)
plt.plot(-60+dx/2+arange(nx)*dx,dpiaLand[:,0,0])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaLand[:,0,1])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaLand[:,0,2],'*')
plt.ylabel('dB')
plt.title('Bright Band')
f1.axes.get_xaxis().set_visible(False)
plt.legend(['CMB_radar','DSRT'])
plt.xlim(-60,60)
plt.subplot(212)
plt.plot(-60+dx/2+arange(nx)*dx,dpiaLand[:,1,0])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaLand[:,1,1])
plt.plot(-60+dx/2+arange(nx)*dx,dpiaLand[:,1,2],'*')
plt.ylabel('dB')
plt.title('No Bright Band')
plt.xlabel('Latitude')
plt.xlim(-60,60)
plt.savefig('landStZonalDPIAGN_SRT.png')

