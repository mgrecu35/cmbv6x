using PyCall
nmu=Ref{Int32}(5)
nmfreq=Ref{Int32}(8)

ccall((:readtablesliang2_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
      nmu,nmfreq)

pushfirst!(PyVector(pyimport("sys")."path"), "")
m=pyimport("readComb")
np=pyimport("numpy")

tb_data=np.loadtxt("tbOut2")
qv,sfcTemp,press,wc,simTb_CMB,pType,envNode,
Nw,psdNodes,binNodes,airTemp,sfcEmiss,sfcBin,cwc=m.readcomb2()
println(size(sfcBin))

tbL=[]
tbL2=[]
tbL3=[]
freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
iFreq=[0,0,1,1,2,3,3,4,4,5,5,6,7]
nfreq=8
tbL1=zeros(0)
tbL2=zeros(0)
kextH3d=zeros(Float32,300,25,88,8)
salbH3d=zeros(Float32,300,25,88,8)
asymH3d=zeros(Float32,300,25,88,8)
for i=1:300
    for j=1:25
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
            for k=11:min(binNodes[i,j,3]+2,88)
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
            for k=binNodes[i,j,3]+3:min(binNodes[i,j,5]+1,88)
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
            kextH3d[i,j,:,:]=kextH
            salbH3d[i,j,:,:]=salbH
            asymH3d[i,j,:,:]=asymH
        end
    end
end
dx=0.25
dx0=5.
idir=1

for i=10:300-10
    for j=12:12
        kextH=zeros(Float32,88,8)
        salbH=zeros(Float32,88,8)
        asymH=zeros(Float32,88,8)
      
        if pType[i,j]>=1
            for k=11:min(binNodes[i,j,5]+1,88)
                i1=Int32(floor(i+((88-k)*dx/dx0+0.5)*idir))
                kextH[k,:]=kextH3d[i1,j,k,:]
                salbH[k,:]=salbH3d[i1,j,k,:]
                asymH[k,:]=asymH3d[i1,j,k,:]
                #println("$i $(i1)")
            end
            #println(binNodes[i,j,:])
            k1=binNodes[i,j,5]+1
            for k=binNodes[i,j,5]+2:sfcBin[i,j]+1
                kextH[k,:]=kextH[k1,:]
                salbH[k,:]=salbH[k1,:]
                asymH[k,:]=asymH[k1,:]
                #println(kextH[k,:])
            end
            #exit(0)
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
                kextInt=np.interp(bins,envNode[i,j,:].+1,kextL)
                tempInt=zeros(88).+273.15
                tempInt[bins]=np.interp(bins,envNode[i,j,:],airTemp[i,j,:])
                z_clw=Array{Float32}(undef,1)
                for k=envNode[i,j,1]+1:envNode[i,j,end]+1
                    if cwc[i,j,k]>0
                        ccall((:gcloud_,"combAlg"),Cvoid,(Ref{Float32},Ref{Float32},
                                                          Ref{Float32},Ref{Float32}),
                              f,Float32(tempInt[k]),Float32(cwc[i,j,k]),
                              z_clw)
                        k1=k-envNode[i,j,1]
                        kextInt[k1]=kextInt[k1]+z_clw[1];
                        #println(z_clw," ",f," ",tempInt[k])
                        #exit(0)
                    end
                end
                kextInt_old=copy(kextInt);

                kextInt=kextInt+kextH[envNode[i,j,1]+1:envNode[i,j,end]+1,iFreq[iEnum]+1]
                salb=copy(kextInt)*0.0;
                asym=copy(kextInt)*0.0;
                fred=(1.0.-kextInt_old./kextInt);
                salb=fred.*copy(salbH[envNode[i,j,1]+1:envNode[i,j,end]+1,iFreq[iEnum]+1]);             
                asym=copy(asymH[envNode[i,j,1]+1:envNode[i,j,end]+1,iFreq[iEnum]+1])
                tLayer=np.interp(bins,envNode[i,j,:],airTemp[i,j,:])
                append!(tLayer,sfcTemp[i,j])
                
                #emis,ebar=pyHB2.emit(f,pol,sfcTemp[i,j],w10[i,j],umu)
                nL=Int32(size(kextInt)[1]);
                umu=cos(53.0/180.0*3.1415);
                umu=Ref{Float32}(umu);
                tb=Array{Float32}(undef,1);
                hght=Float32.((0:nL).*0.25);
                kext=Float32.(kextInt[end:-1:1]);
                salb=Float32.(salb[end:-1:1]);
                asym=Float32.(asym[end:-1:1]);
                fisot=Float32(2.7);
                emis=Float32(sfcEmiss[i,j,iEnum]);
                ebar=Float32(sfcEmiss[i,j,iEnum]);
                
                tLayer=Float32.(tLayer[end:-1:1]);
                sfcTemp1=Float32(sfcTemp[i,j]);
           
                ccall((:radtran2_,"combAlg"), Cvoid,(Ref{Float32},Ref{Int32},Ref{Float32},
                                                       Ref{Float32},Ref{Float32},Ref{Float32},
                                                       Ref{Float32},Ref{Float32},Ref{Float32},
                                                       Ref{Float32},Ref{Float32},Ref{Float32}),
                      umu,nL,tb,sfcTemp1,tLayer,hght,kext,salb,asym,fisot,emis,ebar)
                if iEnum==13
                    if simTb_CMB[i,j,iEnum]>0
                        println(envNode[i,j,:])
                        println(kextL)
                        println(tb[1], " ", simTb_CMB[i,j,iEnum]);
                        append!(tbL1,tb[1]);
                        append!(tbL2,simTb_CMB[i,j,iEnum]);
                        #if simTb_CMB[i,j,iEnum]<tb[1]-30
                        #    println("$i $j $f")
                        #    exit(0)
                        #end
                    end
                end
                iEnum=iEnum+1;
            end
        end
    end
end


tbSim=[]
   
Float32.((1:10)*0.25)
