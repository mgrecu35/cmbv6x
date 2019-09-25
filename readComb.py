from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import *

def readcomb2():
    fname='2B.GPM.DPRGMI.20161001.026079.HDF5'
    fpath='./'
    fh=Dataset(fpath+fname)
    
    Lat=fh['NS']['Latitude'][:,:]
    Lon=fh['NS']['Longitude'][:,:]
    
    a=[range(300)]
    m = Basemap(width=400000,height=400000,
                rsphere=(6378137.00,6356752.3142),\
                    resolution='l',area_thresh=1000.,projection='lcc',\
                    lat_1=-3.,lat_2=1,lat_0=-1,lon_0=97.)
                
    m.drawcoastlines()
    x,y=m(Lon[a[0],:],Lat[a[0],:])
    
    sfcRain=fh['NS/surfPrecipTotRate'][a[0],:]
    #plt.contourf(x,y,sfcRain,levels=[0.1,0.2,0.4,1,2,4,8,16,32,64],vmin=0,cmap='jet')
    qv=fh['NS']['vaporDensity'][a[0],:,:]
    qv_sfc=fh['NS/surfaceVaporDensity'][a[0],:]
    simTb_CMB=fh['NS/simulatedBrightTemp'][a[0],:,:]
    envNode=fh['NS/envParamNode'][a[0],:,:]
    sfcTemp=fh['NS/skinTemperature'][a[0],:]
    airTemp=fh['NS/airTemperature'][a[0],:,:]
    press=fh['NS/airPressure'][a[0],:,:]
    sfcEmiss=fh['NS/surfEmissivity'][a[0],:,:]
    w10=fh['NS/tenMeterWindSpeed'][a[0],:]
    pType=(fh['NS/Input/precipitationType'][a[0],:]/1e7).astype(int)
    binNodes=fh['NS/phaseBinNodes'][a[0],:,:]
    Nw=fh['NS/precipTotPSDparamLow'][a[0],:,:]-log10(8e6)
    psdNodes=fh['NS/PSDparamLowNode'][a[0],:,:]
    wc=fh['NS/precipTotWaterCont'][a[0],:,:]
    return qv,sfcTemp,press,wc,simTb_CMB,pType,envNode,Nw,psdNodes,binNodes,airTemp,sfcEmiss
