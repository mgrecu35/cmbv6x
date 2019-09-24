import matplotlib
from netCDF4 import Dataset
gui='tkagg'
matplotlib.use(gui,warn=False, force=True)
import matplotlib.pyplot as plt
import matplotlib.colors as col
from numpy import *
def pcmesh(z,zmin,zmax):
    print(z.shape)
    print(type(z))
    print(z)
    print(zmin,zmax)
    plt.figure()
    z=array(z)
    print(type(z))
    plt.pcolormesh(z,vmin=zmin,vmax=zmax,cmap='jet')
    #plt.show()

def pshow():
    plt.show()
z=random.randn(50,50)
pcmesh(z,-1,1)
