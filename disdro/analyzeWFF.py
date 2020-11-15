
from numpy import *
from bhmief import *
import pytmatrix.refractive
wl=[pytmatrix.refractive.wl_Ku,\
    pytmatrix.refractive.wl_Ka,pytmatrix.refractive.wl_W]
#print(wl) # units are mm
#stop

mw1= pytmatrix.refractive.m_w_10C[wl[0]]
mw2= pytmatrix.refractive.m_w_10C[wl[1]]
mw3= pytmatrix.refractive.m_w_10C[wl[2]]

sback13w=[]
sback37w=[]
sback95w=[]
nang=20

dbinsM=(0.025+arange(300)*0.05)
massW=(0.1*dbinsM)**3/6.*pi
dMax=10.*(massW/0.0061)**(1/2.05)

mi1= pytmatrix.refractive.mi(wl[0],0.92)
mi2= pytmatrix.refractive.mi(wl[1],0.92)
mi3= pytmatrix.refractive.mi(wl[2],0.92)
m=array([mi1,mi2,mi3])

K2=0.93
kextw_13=[]
kextw_37=[]
kextw_95=[]
vdop_13=[]
vdop_37=[]
vdop_95=[]
dbinsM=arange(50)*0.2+0.1
vT=3.78*dbinsM**0.67

for i,d in enumerate(dbinsM):
    x=d*pi/wl[0]
    s1,s2,qext13w,qsca13w,qback13w,gsca13w=bhmie(x,mw1,nang)
    sback13w.append(qback13w*d**2*pi/4.*wl[0]**4/pi**5) ## mm^6
    kextw_13.append(qext13w*d**2*pi/4.) ##mm^2
    x=d*pi/wl[1]
    s1,s2,qext37w,qsca37w,qback37w,gsca37w=bhmie(x,mw2,nang)
    sback37w.append(qback37w*(d)**2*pi/4.*wl[1]**4/pi**5)
    kextw_37.append(qext37w*d**2*pi/4.) ##mm^2
    x=d*pi/wl[2]
    s1,s2,qext95w,qsca95w,qback95w,gsca95w=bhmie(x,mw3,nang)
    sback95w.append(qback95w*(d)**2*pi/4.*wl[2]**4/pi**5)
    kextw_95.append(qext95w*d**2*pi/4.) ##mm^2

sback13w=array(sback13w)
sback37w=array(sback37w)
sback95w=array(sback95w)
kextw_13=array(kextw_13)
kextw_37=array(kextw_37)
kextw_95=array(kextw_95)

massW=(0.1*dbinsM)**3/6.*pi
def procWFF(lines):
    Z=[]
    att=[]
    dm=[]
    lwc=[]
    y=[]
    NwL=[]
    vDop=[]
    for l in lines:
        ls=l.split()
        Nd=array([float(nc) for nc in ls[4:]])*0.2
        zKu=log10(sum(Nd*sback13w/K2)+1e-9)*10.
        zKa=log10(sum(Nd*sback37w/K2)+1e-9)*10.
        zW=log10(sum(Nd*sback95w/K2)+1e-9)*10.
        vDopKu=sum(Nd*sback13w*vT)/sum(Nd*sback13w)
        vDopKa=sum(Nd*sback37w*vT)/sum(Nd*sback37w)
        #dm=sum(Nd*dbinsM**4)/sum(Nd*dbinsM**3)
        kextKu=sum(Nd*kextw_13)*1e-3
        kextKa=sum(Nd*kextw_37)*1e-3
        kextW=sum(Nd*kextw_95)*1e-3
        M=sum(Nd*massW)
        lwc.append(M)
        dm1=sum(Nd*massW*dbinsM)/M
        dm.append(dm1)
        Nw=4**4/pi*M/(0.1*dm1)**4/1e6
        if int(ls[1])>120 and int(ls[1])<300 :
            print(ls[0])
            Z.append([zKu,zKa,zW])
            att.append([kextKu,kextKa,kextW])
            y.append([M,dm1,Nw])
            vDop.append([vDopKu,vDopKa,dm1])
    return array(Z), array(att), array(y),array(vDop)
from getZ import *
f='wff_2dvd_sn36_raindsd_ter.txt'



#f='iphex/iphex_2dvd_sn25_20140611_N351335.71_W820323.41_rainDSD.txt'

lines=[]

fs=['iphex/iphex_2dvd_sn25_20140530_N351335.71_W820323.41_rainDSD.txt',
    'iphex/iphex_2dvd_sn25_20140611_N351335.71_W820323.41_rainDSD.txt',
    'iphex/iphex_2dvd_sn35_20140515_N351734.29_W821014.52_rainDSD.txt',
    'wff_2dvd_sn25_raindsd_ter.txt',
    'wff_2dvd_sn35_raindsd_ter.txt',
    'wff_2dvd_sn36_raindsd_ter.txt',
    'wff_2dvd_sn37_raindsd_ter.txt',
    'wff_2dvd_sn38_raindsd_ter.txt',
    'wff_2dvd_sn70_raindsd_ter.txt']
for f in fs:
    lines.extend(open(f).readlines())
zObs1,atten1,y1,vdop=procWFF(lines)



import matplotlib.pyplot as plt
Nw1=y1[:,-1]
dm=y1[:,1]
a=nonzero(y1[:,0]<1)
b=nonzero(y1[:,0][a]>0.05)
plt.scatter((dm[a][b]),log10(Nw1[a][b]))
zd=[]
for i in range(15*4):
    di1=0.2+i*0.125/2
    di2=0.2+(i+1)*0.125/2
    c=nonzero((dm[a][b]-di1)*(dm[a][b]-di2)<0)
    if len(c[0]>0):
        print(zObs1[a[0][b][c],1].mean(),dm[a][b][c].mean())
        zd.append([zObs1[a[0][b][c],0].mean(),dm[a][b][c].mean()])

plt.figure()
plt.scatter(array(zd)[:,0],array(zd)[:,1])
polycoeff=polyfit(array(zd)[:,0],log(array(zd)[:,1]),2)
plt.figure()
a=nonzero(zObs1[:,0]>-10)
plt.scatter(dm[a[0]],log10(4.343*atten1[a[0],1]/10**(0.1*zObs1[a[0],0])))
import pickle
x1,y1=pickle.load(open("../cmbT.pklz","rb"))
plt.plot(1.00*x1,y1,color='red')
plt.legend(['CMB','Disdrometer'])
plt.xlim(0.1,4)
plt.xlabel("D$_m$(mm)")
plt.ylabel("log10(k$_{Ka}$/Z$_{Ku}$)")
plt.title("Combined Algorithm & Disdrometer Tables")
plt.savefig("kOverZ_dm.png")


stop
nwm=Nw1.mean()
for i in range(1):
    x=zObs1[:,i]/10.-log10(Nw1/nwm)
    y=log10(4.4343*atten1[:,i]/(Nw1/nwm))
    res=polyfit(x,y,1)
    print(res)
    plt.figure()
    plt.scatter(zObs1[:,i]-log10(Nw1/nwm)*10,log10(y1[:,0]/y1[:,-1]*nwm))

plt.figure()
import matplotlib
matplotlib.rcParams.update({'font.size': 13})
rratio=10.**y/10**(res[0]*x)/10**(res[1])
log_att=polyval(res,x)
pratio=10.**log_att/10**(res[0]*x)/10**(res[1])
plt.scatter(10*x,rratio)
plt.scatter(10*x,pratio)
plt.xlabel('Z-10*log10(N$_w$/N$_{wmean}$) (dBZ)')
plt.ylabel('k/log10(N$_w$)/$\\alpha$/(Z N$_{wmean}$/N$_w$) $^\\beta$')
plt.legend(['True','Power law'])
plt.savefig('kZ.png')
stop
    
