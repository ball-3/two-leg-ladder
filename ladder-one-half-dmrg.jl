using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(6969)

function main(length, jl, jr)
	
	N = length*2
	Jl = jl		#inter dimer interaction
	Jr = jr		#intra dimer interaction
			#we are testing limit Jr >> Jl

	sites = siteinds("S=1/2", N)
	
	os = OpSum()
	for j in 1:(N-2)
		#inter dimer
		os += Jl, "Sz", j, "Sz", j+2
		os += Jl*0.5, "S+", j, "S-", j+2
		os += Jl*0.5, "S-", j, "S+", j+2
	end
	
	for j in 1:2:(N-1)	
		#intra dimer
		os += Jr, "Sz", j, "Sz", j+1
		os += Jr*0.5, "S+", j, "S-", j+1
		os += Jr*0.5, "S-", j, "S+", j+1
	end
	
	H = MPO(os, sites)
	
	psi0 = randomMPS(sites; linkdims=10)
	
	nsweeps = 20#changed from 10
	maxdim = [200]#increased
	
	cutoff = [1e-12]#changed from 1e-6
	
	energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
	
	#@printf("Final energy = %.12f, n = %d, Jl = %f, Jr = %f\n", energy, N, Jl, Jr)
	
	return [energy, psi]
end

#main(4,1,1)

