using PyCall

pushfirst!(PyVector(pyimport("sys")."path"), "")
m=pyimport("readComb")
qv,sfcTemp,press,wc=m.readcomb2()