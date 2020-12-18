iread=1
from numpy import *
import matplotlib.pyplot as plt
import pickle
import glob

fs=sorted(glob.glob("retrs/stratRetrievalsAug2018StE*"))
dx=5
nx=int(120/dx)
ny=int(360/dx)
countsLand=zeros((nx,2),float)
countsOcean=zeros((nx,2),float)
rainLand=zeros((nx,2,4),float)
rainOcean=zeros((nx,2,4),float)
dpiaLand=zeros((nx,2,3),float)
dpiaOcean=zeros((nx,2,3),float)
dmLand=zeros((nx,2,3),float)
dmOcean=zeros((nx,2,3),float)

counts2d=zeros((nx,ny,2),float)
rain2d=zeros((nx,ny,2,4),float)
dpia2d=zeros((nx,ny,2,3),float)
dm2d=zeros((nx,ny,2,3),float)


rs=[]
rsL=[]
zsfcL=[]
zsfcO=[]
dpiaO=[]
dpiaL=[]
for f1 in fs[:]:
    [r1L,r1Lun]=pickle.load(open(f1,'rb'))

    for r1 in r1L:
        rate,dprrate,dpia,dpiaSRT,sfcType,binBB,lat,lon,cmbSfcRain,\
            rate_noad,dpia_noad,zKa,zkasim,zkasim_noad,\
            bzd,binBB,bcf,rrate1df,zkasimF,dpiaF,dmDPR,dmCMB,dm_New=r1
        
        if type(dmDPR)!=float32:
            continue
        i0=int((lat+60)/dx)
        i0=max(0,i0)
        i0=min(i0,nx-1)
        j0=int((lon+180)/dx)
        j0=max(0,j0)
        j0=min(j0,ny-1)
        if binBB>0:
            ibb=0
        else:
            ibb=1
        if sfcType==0:
            rs.append([rate,dprrate,cmbSfcRain,rrate1df])
            dpiaO.append([dpia,dpiaSRT,dpiaF])
            if zKa[bcf]>10 and binBB>0:
                zsfcO.append([zKa[bcf],zkasim[bcf],zkasim_noad[bcf],zkasimF])
            countsOcean[i0,ibb]+=1
            rainOcean[i0,ibb,0]+=rate
            rainOcean[i0,ibb,3]+=rrate1df
            rainOcean[i0,ibb,1]+=dprrate
            dpiaOcean[i0,ibb,0]+=dpia
            dpiaOcean[i0,ibb,2]+=dpiaF
            dpiaOcean[i0,ibb,1]+=dpiaSRT
            dmOcean[i0,ibb,0]+=dmDPR
            dmOcean[i0,ibb,2]+=dmCMB
            dmOcean[i0,ibb,1]+=dm_New
            rainOcean[i0,ibb,2]+=cmbSfcRain
        else:
            rsL.append([rate,dprrate,cmbSfcRain,rrate1df])
            dpiaL.append([dpia,dpiaSRT,dpiaF])
            if zKa[bcf]>10 and binBB>0:
                zsfcL.append([zKa[bcf],zkasim[bcf],zkasim_noad[bcf],zkasimF])
            countsLand[i0,ibb]+=1
            rainLand[i0,ibb,0]+=rate
            rainLand[i0,ibb,3]+=rrate1df
            rainLand[i0,ibb,1]+=dprrate
            dpiaLand[i0,ibb,0]+=dpia
            dpiaLand[i0,ibb,2]+=dpiaF
            dpiaLand[i0,ibb,1]+=dpiaSRT
            rainLand[i0,ibb,2]+=cmbSfcRain
            dmLand[i0,ibb,0]+=dmDPR
            dmLand[i0,ibb,2]+=dmCMB
            dmLand[i0,ibb,1]+=dm_New
        rain2d[i0,j0,ibb,3]+=rate
        rain2d[i0,j0,ibb,1]+=rrate1df
        rain2d[i0,j0,ibb,0]+=dprrate
        rain2d[i0,j0,ibb,2]+=cmbSfcRain
        dpia2d[i0,j0,ibb,0]+=dpia
        dpia2d[i0,j0,ibb,2]+=dpiaF
        dpia2d[i0,j0,ibb,1]+=dpiaSRT
        dm2d[i0,j0,ibb,0]+=dmDPR
        dm2d[i0,j0,ibb,2]+=dmCMB
        dm2d[i0,j0,ibb,1]+=dm_New
        counts2d[i0,j0,ibb]+=1
        
a=(countsLand>0)
rmax=6
rainLand[:,:,0][a]/=countsLand[a]
rainLand[:,:,1][a]/=countsLand[a]
dpiaLand[:,:,0][a]/=countsLand[a]
dpiaLand[:,:,1][a]/=countsLand[a]
dpiaLand[:,:,2][a]/=countsLand[a]
dmLand[:,:,0][a]/=countsLand[a]
dmLand[:,:,1][a]/=countsLand[a]
dmLand[:,:,2][a]/=countsLand[a]
rainLand[:,:,2][a]/=countsLand[a]
rainLand[:,:,3][a]/=countsLand[a]
a=(countsOcean>0)
rainOcean[:,:,0][a]/=countsOcean[a]
rainOcean[:,:,1][a]/=countsOcean[a]
dpiaOcean[:,:,0][a]/=countsOcean[a]
dpiaOcean[:,:,1][a]/=countsOcean[a]
dpiaOcean[:,:,2][a]/=countsOcean[a]
dmOcean[:,:,0][a]/=countsOcean[a]
dmOcean[:,:,1][a]/=countsOcean[a]
dmOcean[:,:,2][a]/=countsOcean[a]
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



plt.figure()
plt.suptitle('Ocean')
f1=plt.subplot(211)
plt.plot(-60+dx/2+arange(nx)*dx,dmOcean[:,0,0])
plt.plot(-60+dx/2+arange(nx)*dx,dmOcean[:,0,1])
plt.plot(-60+dx/2+arange(nx)*dx,dmOcean[:,0,2])
plt.ylabel('mm')
plt.title('Bright Band')
f1.axes.get_xaxis().set_visible(False)
plt.legend(['DPR','CMB_V6','CMB_radar'])
plt.xlim(-60,60)
plt.subplot(212)
plt.plot(-60+dx/2+arange(nx)*dx,dmOcean[:,1,0])
plt.plot(-60+dx/2+arange(nx)*dx,dmOcean[:,1,1])
plt.plot(-60+dx/2+arange(nx)*dx,dmOcean[:,1,2])
plt.ylabel('mm')
plt.title('No Bright Band')
plt.xlabel('Latitude')
plt.xlim(-60,60)
plt.savefig('dmOceanStZonal_SRT.png')


plt.figure()
plt.suptitle('Land')
f1=plt.subplot(211)
plt.plot(-60+dx/2+arange(nx)*dx,dmLand[:,0,0])
plt.plot(-60+dx/2+arange(nx)*dx,dmLand[:,0,1])
plt.plot(-60+dx/2+arange(nx)*dx,dmLand[:,0,2])
plt.ylabel('mm')
plt.title('Bright Band')
f1.axes.get_xaxis().set_visible(False)
plt.legend(['DPR','CMB_V6','CMB_radar'])
plt.xlim(-60,60)
plt.ylim(0.8,1.6)
plt.subplot(212)
plt.plot(-60+dx/2+arange(nx)*dx,dmLand[:,1,0])
plt.plot(-60+dx/2+arange(nx)*dx,dmLand[:,1,1])
plt.plot(-60+dx/2+arange(nx)*dx,dmLand[:,1,2])
plt.ylabel('mm')
plt.title('No Bright Band')
plt.xlabel('Latitude')
plt.xlim(-60,60)
plt.ylim(0.8,1.6)
plt.savefig('dmLandStZonal_SRT.png')


a=nonzero(counts2d>0)
b=nonzero(counts2d==0)


for i in range(3):
    dm2d[:,:,:,i][a]/=counts2d[a]

for i in range(4):
    rain2d[:,:,:,i][a]/=counts2d[a]
    

from mpl_toolkits.basemap import Basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=60,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
xx=-180+dx*arange(ny+1)+dx/2
yy=-60+dx*arange(nx+1)+dx/2
x2,y2=meshgrid(xx,yy)
x2m,y2m=m(x2,y2)




plt.figure(figsize=(11,8))

dm2dL=[]
for i in range(3):
    dm2dL.append(ma.array(dm2d[:,:,0,i],mask=counts2d[:,:,0]==0))
    
plt.suptitle('Stratiform with BB')
f1=plt.subplot(311)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,dm2dL[0],\
               cmap='jet',vmin=0.8,vmax=1.6)
plt.title('DPR')
plt.colorbar()
f1=plt.subplot(312)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,dm2dL[2],\
               cmap='jet',vmin=0.8,vmax=1.6)
plt.title('CMB_V6')
plt.colorbar()
f1=plt.subplot(313)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,dm2dL[1],\
               cmap='jet',vmin=0.8,vmax=1.6)
plt.title('CMB_Rad')
plt.colorbar()
plt.savefig('dmAug2018StWithBB.png')
                

plt.figure(figsize=(11,8))
plt.suptitle('Stratiform without BB')
dm2dL=[]
for i in range(3):
    dm2dL.append(ma.array(dm2d[:,:,1,i],mask=counts2d[:,:,1]==0))

f1=plt.subplot(311)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,dm2dL[0],\
               cmap='jet',vmin=0.8,vmax=1.6)
plt.title('DPR')
plt.colorbar()
f1=plt.subplot(312)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,dm2dL[2],\
               cmap='jet',vmin=0.8,vmax=1.6)
plt.title('CMB_V6')
plt.colorbar()
f1=plt.subplot(313)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,dm2dL[1],\
               cmap='jet',vmin=0.8,vmax=1.6)
plt.title('CMB_Rad')
plt.colorbar()
plt.savefig('dmAug2018StWithOutBB.png')
                
                


plt.figure(figsize=(11,8))

rain2dL=[]
for i in range(3):
    rain2dL.append(ma.array(rain2d[:,:,0,i],mask=counts2d[:,:,0]==0))
    
plt.suptitle('Stratiform with BB')
f1=plt.subplot(311)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,rain2dL[0],\
               cmap='jet',vmin=1,vmax=8)
plt.title('DPR')
plt.colorbar()
f1=plt.subplot(312)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,rain2dL[2],\
               cmap='jet',vmin=1,vmax=8)
plt.title('CMB_V6')
plt.colorbar()
f1=plt.subplot(313)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,rain2dL[1],\
               cmap='jet',vmin=1,vmax=8)
plt.title('CMB_Rad')
plt.colorbar()
plt.savefig('rainAug2018StWithBB.png')
                

plt.figure(figsize=(11,8))
plt.suptitle('Stratiform without BB')
rain2dL=[]
for i in range(3):
    rain2dL.append(ma.array(rain2d[:,:,1,i],mask=counts2d[:,:,1]==0))

f1=plt.subplot(311)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,rain2dL[0],\
               cmap='jet',vmin=1,vmax=6)
plt.title('DPR')
plt.colorbar()
f1=plt.subplot(312)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,rain2dL[2],\
               cmap='jet',vmin=1,vmax=6)
plt.title('CMB_V6')
plt.colorbar()
f1=plt.subplot(313)
m.drawcoastlines()
plt.pcolormesh(x2m,y2m,rain2dL[1],\
               cmap='jet',vmin=1,vmax=6)
plt.title('CMB_Rad')
plt.colorbar()
plt.savefig('rainAug2018StWithOutBB.png')
                
                
