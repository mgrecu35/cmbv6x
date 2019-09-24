using PyCall
push!(PyVector(pyimport("sys")["path"]), ".")
p=pyimport("pyplot")