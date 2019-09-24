it=Array{Int32}(undef,1)
ccall((:testr_,"combAlg"),Cvoid,(Ref{Int32},),it)
println(it[1])
#exit(0)
fname="paramFile"
ifs=Ref{Int32}(1)
ifs=Array{Int32,1}(undef,1)
ifs[1]=1
ccall((:mainj,"combAlg"),Int32,(Int32,Cstring,Ref{Int32}),1,fname,ifs)

ny=49
nx=300
nz=88
zku=Array{Float32}(undef,nx*ny*nz)	
zka=Array{Float32}(undef,nx*ny*nz)	
rainType=Array{Int32}(undef,nx*ny)	

#using PyPlot
using PyCall
push!(PyVector(pyimport("sys")["path"]), ".")
p=pyimport("pyplot")
#pygui(true)
#zku3d=reshape(zku,nz,ny,nx)


tbRgrid=Array{Float32}(undef,14*49*9300)
oe_tbs=Array{Float32}(undef,300*49*13*4)
scattProp=Array{Float32}(undef,13*300*49*88*3)
dprrain=Array{Float32}(undef,49*300) 
println("ifs=$(ifs[1])")
#exit(0)
for i1=0:6
	i=Int32(i1)
	idir=Int32(1)
	if(ifs[1]!=1) 
		ccall((:do_chunk_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32},Ref{Int32}),i,1,idir)
	else
		ccall((:do_chunkx_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32},Ref{Int32}),i,1,idir)
	end
	if i>=4
		println(idir)
		#exit(0);
		idir=1
		ccall((:radarretsub2_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32},#
		Ref{Int32},Ref{Float32},Ref{Float32},Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32}),
		5,8,i*300,tbRgrid,dprrain,i,100,1,idir)
		ccall((:radarretsub3_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32},#
		Ref{Int32},Ref{Float32},Ref{Float32},Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32}),
		5,8,i*300,tbRgrid,dprrain,i,100,1,idir)
		if(ifs[1]==1) 
	            ccall((:radarretsub4_fs_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32},#
		    Ref{Int32},Ref{Float32},Ref{Float32},Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32}),
	            5,8,i*300,tbRgrid,dprrain,i,100,1,idir)
		else
                    ccall((:radarretsub4_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32},#
		    Ref{Int32},Ref{Float32},Ref{Float32},Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32}),
	            5,8,i*300,tbRgrid,dprrain,i,100,1,idir)		    
		end
		nscans=Int32(300)
		ccall((:get_rain_type_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),rainType,it)
		rainType2d=reshape(rainType,ny,300)
		nscans=it[1]
		println(nscans)
		
		
		ccall((:get_oetbs_,"combAlg"),Cvoid,(Ref{Float32},Ref{Float32}),oe_tbs,scattProp)
		ccall((:dealloc_struct_,"combAlg"),Cvoid,(Ref{Int32},),i)
		
	end
	ccall((:getzs_,"combAlg"),Cvoid,(Ref{Float32},Ref{Float32},Ref{Int32},Ref{Int32},Ref{Int32}),	zku,zka,nx,ny,nz)	 	
	ccall((:dealloc_chunk_,"combAlg"),Cvoid,(Ref{Int32},),i)
	
end

ccall((:closefiles_,"combAlg"),Cvoid,(Ref{Int32},),1)
oe_tbs=reshape(oe_tbs,300,49,13,4)
scattProp=reshape(scattProp,13,300,49,88,3)
xr=pyimport("xarray")
scxr=xr.DataArray(scattProp);
scdset=xr.Dataset(Dict("scattProp"=>scxr))
scdset.to_netcdf("scattProp.nc")

tbRgrid=reshape(tbRgrid,14,49,9300)[:,:,1801:2100]
zku3d=reshape(zku,nz,ny,nx)
zka3d=reshape(zka,nz,ny,nx)
zx=copy(zku3d[end:-1:1,24,:])
p.pcmesh(copy(tbRgrid[7,12:37,:]'),100,270)
#p.pcmesh(zx,0.0,50.0)
p.pshow()
using Statistics
a=findall(oe_tbs[:,:,7,2].>0)

#for i1=1:4 
#	println(cor(tbRgrid[13,13:37,:]'[a],oe_tbs[:,:,13,i1][a]))
#end	
nxs=275
zx1=copy(zku3d[end:-1:1,:,nxs]);
#p.pcmesh(zx1,0,50)
zx2=copy(zka3d[end:-1:1,:,nxs]);
#p.pcmesh(zx2,0,50)
zxk=copy(zka3d[end:-1:1,24,:])
#p.pcmesh(zxk[end:-1,24,200:300],0,50)
