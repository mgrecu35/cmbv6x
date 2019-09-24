f=netcdf.Dataset("2B.GPM.DPRGMI.20161001.026079.HDF5")
wc=get(get(f,"MS"),"precipTotWaterCont");