import pickle
import matplotlib.pyplot as plt
from numpy import *
dataL=pickle.load(open('sp_noG_Data.pklz','rb'))

latL=[]
lonL=[]
counts=zeros((120,3),float)
sfcrain=zeros((120,3),float)
dprsfcrain=zeros((120,3),float)
dsrtPIAL=zeros((120,2),float)
for d1 in dataL:
    rrate1d,dprsfcRate,bst,bcf,lon,lat,\
        reliabFlag,dsrtPIA,piaka,piaku,zKu1,pType,sfcType=d1
    latL.append(lat)
    lonL.append(lon)
    #if reliabFlag==0:
    #    print(pType)
    #if sfcType==0:
    #    continue
    if lat>-60 and lat<60:
        if pType==1 and reliabFlag==1:
            i0=int(lat-60)
            sfcrain[i0,0]+=rrate1d[bcf]
            dprsfcrain[i0,0]+=dprsfcRate
            counts[i0,0]+=1
            dsrtPIAL[i0,0]+=dsrtPIA
            dsrtPIAL[i0,1]+=piaka-piaku
        if pType==1 and reliabFlag!=1:
            i0=int(lat-60)
            sfcrain[i0,1]+=rrate1d[bcf]
            dprsfcrain[i0,1]+=dprsfcRate
            counts[i0,1]+=1
        if pType==2:
            i0=int(lat-60)
            sfcrain[i0,2]+=rrate1d[bcf]
            dprsfcrain[i0,2]+=dprsfcRate
            counts[i0,2]+=1

a=nonzero(counts>0)
sfcrain[a]/=counts[a]
dprsfcrain[a]/=counts[a]
plt.figure(figsize=(8,8))
plt.subplot(311)
plt.plot(-59.5+arange(120),sfcrain[:,0])
plt.plot(-59.5+arange(120),dprsfcrain[:,0])
plt.ylim(0,6)
plt.xlim(-60,60)
plt.subplot(312)
plt.plot(-59.5+arange(120),sfcrain[:,1])
plt.plot(-59.5+arange(120),dprsfcrain[:,1])
plt.ylim(0,6)
plt.xlim(-60,60)
plt.subplot(313)
plt.plot(-59.5+arange(120),sfcrain[:,2])
plt.plot(-59.5+arange(120),dprsfcrain[:,2])
plt.ylim(0,30)
plt.xlim(-60,60)

plt.figure()
plt.plot(dsrtPIAL[:,1]/counts[:,0])
plt.plot(dsrtPIAL[:,0]/counts[:,0])
