using PyCall
nmu=Ref{Int32}(5)
nmfreq=Ref{Int32}(8)

ccall((:readtablesliang2_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
      nmu,nmfreq)

pushfirst!(PyVector(pyimport("sys")."path"), "")
m=pyimport("readComb")
np=pyimport("numpy")
qv,sfcTemp,press,wc,simTb_CMB,pType,envNode,
Nw,psdNodes,binNodes,airTemp,sfcEmiss=m.readcomb2()

tbL=[]
tbL2=[]
tbL3=[]
freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
iFreq=[0,0,1,1,2,3,3,4,4,5,5,6,6,7,7]
nfreq=8
for i=1:300
    for j=1:49
        kextH=zeros(Float32,88,8)
        salbH=zeros(Float32,88,8)
        asymH=zeros(Float32,88,8)
        Nwb=zeros(88)
        if pType[i,j]>=1
            xnw=1:88;
            xnw=(psdNodes[i,j,1]+1:psdNodes[i,j,end]+1);
            Nwb[psdNodes[i,j,1]+1:psdNodes[i,j,end]+1]=
            np.interp(xnw,
                      psdNodes[i,j,:],Nw[i,j,:,1])
            for k=binNodes[i,j,1]+1:min(binNodes[i,j,3]+1,88)
                n0w=Ref{Float32}(Nwb[k])
                iwc=Ref{Float32}(wc[i,j,k])
                kext=Array{Float32}(undef,8)
                salb=Array{Float32}(undef,8)
                asym=Array{Float32}(undef,8)
                ccall((:getsnowp_2_,"combAlg"),Cvoid,(Ref{Float32},Ref{Float32},
                                                     Ref{Float32},Ref{Float32},
                                                     Ref{Float32},Ref{Int32}),
                     n0w,iwc,kext,salb,asym,nmfreq)
                kextH[k,:]=kext
                salbH[k,:]=salb
                asymH[k,:]=asym
            end
            for k=binNodes[i,j,3]+2:min(binNodes[i,j,5]+1,88)
                n0w=Ref{Float32}(Nwb[k])
                iwc=Ref{Float32}(wc[i,j,k])
                kext=Array{Float32}(undef,8)
                salb=Array{Float32}(undef,8)
                asym=Array{Float32}(undef,8)
                ccall((:getrainp_2_,"combAlg"),Cvoid,(Ref{Float32},Ref{Float32},
                                                     Ref{Float32},Ref{Float32},
                                                     Ref{Float32},Ref{Int32}),
                     n0w,iwc,kext,salb,asym,nmfreq)
                kextH[k,:]=kext
                salbH[k,:]=salb
                asymH[k,:]=asym
            end
            iEnum=1
            for (pol,f) in zip(npol,freqs)
                kextL=zeros(Float32,0)
                one=Ref{Int32}(1)
                absair=Array{Float32}(undef,1)
                abswv=Array{Float32}(undef,1)
                for (q,tk,pa) in zip(qv[i,j,:],airTemp[i,j,:],press[i,j,:])
                    ccall((:gasabsr98_,"combAlg"),Cvoid,
                          (Ref{Float32},Ref{Float32},Ref{Float32},
                           Ref{Float32},Ref{Float32},Ref{Float32},Ref{Int32}),
                          f,tk,q*1e-3,pa*1e2,absair,abswv,one)
                    append!(kextL,absair[1]+abswv[1])
                end
                bins=1:2;
                bins=envNode[i,j,1]+1:envNode[i,j,end]+1;
                kextInt=np.interp(bins,envNode[i,j,:],kextL)
                kextInt_old=copy(kextInt);
                kextInt=kextInt+kextH[envNode[i,j,1]+1:envNode[i,j,end]+1,iFreq[iEnum]+1]
                salb=copy(kextInt)*0.0;
                asym=copy(kextInt)*0.0;
                fred=(1.0.-kextInt_old./kextInt);
                salb=fred.*copy(salbH[envNode[i,j,1]+1:envNode[i,j,end]+1,iFreq[iEnum]+1]);             
                asym=copy(asymH[envNode[i,j,1]:envNode[i,j,end]+1,iFreq[iEnum]+1])
                tLayer=np.interp(bins,envNode[i,j,:],airTemp[i,j,:])
                append!(tLayer,sfcTemp[i,j])
                iEnum=iEnum+1
                #emis,ebar=pyHB2.emit(f,pol,sfcTemp[i,j],w10[i,j],umu)
                nL=size(kextInt);
            end
        end
    end
end


tbSim=[]
   
#for k in range(binNodes[i,j,0],\
#                           min(binNodes[i,j,3]+1,binNodes[i,j,4])):
#                kext,salb,asym=pyHB2.getsnowp(0.5+Nwb[k-psdNodes[i,j,0]],\
#                                            pRate[i,j,k],nfreq)
#                kextH[k,:]=kext
#                salbH[k,:]=salb
#                asymH[k,:]=asym
#            for k in range(binNodes[i,j,3]+1,binNodes[i,j,4]):
#                kext,salb,asym=pyHB2.getrainp(Nwb[k-psdNodes[i,j,0]],\
#                                            pRate[i,j,k],nfreq)