import matplotlib.pyplot as plt
import matplotlib
from netCDF4 import Dataset
import numpy as np
from numpy import *

f=Dataset("scattProp.nc")
scattProp=f['scattProp'][:,:,:,:,:]
x1=1800
scattPropm=np.ma.array(scattProp,mask=scattProp<0.001)
plt.subplot(311)
fig=plt.pcolormesh(arange(300)+1800,arange(88)*0.25,scattPropm[5,:,5,::-1,0].T,cmap='jet',vmin=0.1, vmax=3.0, \
                       norm=matplotlib.colors.LogNorm())
plt.xlim(x1+150,x1+210)
fig.axes.get_xaxis().set_visible(False)
plt.title("Extinction")
plt.ylim(0,15)
plt.colorbar()
plt.subplot(312)
fig=plt.pcolormesh(arange(300)+1800,arange(88)*0.25,scattPropm[5,:,5,::-1,1].T,cmap='jet',vmin=0.001)
plt.xlim(x1+150,x1+210)
fig.axes.get_xaxis().set_visible(False)
plt.ylim(0,15)
plt.title("Scattering Albedo")
plt.colorbar()
plt.subplot(313)

fig=plt.pcolormesh(arange(300)+1800,arange(88)*0.25,scattPropm[5,:,5,::-1,2].T,cmap='jet',vmin=0)
plt.xlim(x1+150,x1+210)
plt.ylim(0,15)
plt.title("Asymmetry factor")
plt.colorbar()
plt.savefig("scattProp36.6.png")
