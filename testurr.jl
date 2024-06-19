using ITensors, ITensorMPS

function hi(Jl,Jr,N)
	
	sites = siteinds("S=1/2", N)
	os = OpSum()
	
	for j in 1:(N-2)
		#inter dimer
		os += Jl, "Sz", j, 1, j+2
		os += Jl*0.5, "S-", j, 1, j+2
		os += Jl*0.5, "S-", j, "S+", j+2
	end
	
	for j in 1:2:(N-1)	
		#intra dimer
		#os += "Sz"*"Sz", j, j+1
		os += Jr, "Sz",j,"I",j+1
		os += Jr*0.5, "S+", j, "S-", j+1
		os += Jr*0.5, "S-", j, "S+", j+1
	end
	
	display(os)
		
	op = MPO(os, sites)
	display(op)
end

hi(1,1,2)
