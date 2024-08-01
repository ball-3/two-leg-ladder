using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(9999)

function main2(N,jl)

	sites = siteinds("S=1/2", N)
	
	os = OpSum()
		for j in 1:(N-1)
			os += jl,"Sz", j, "Sz", j+1
			os += jl*0.5, "S+", j, "S-", j+1
			os += jl*0.5, "S-", j, "S+", j+1
		end
	
	operator = MPO(os, sites)
	
	psi0 = randomMPS(sites; linkdims=10)
	
	nsweeps = 10
	maxdim = [10,20,100,100,200]
	
	cutoff = [1e-6]
	
	energy, psi = dmrg(operator, psi0; nsweeps, maxdim, cutoff)
	
	return [energy, psi]
end
