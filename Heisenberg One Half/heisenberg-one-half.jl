using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(9999)

function main(N)

	sites = siteinds("S=1/2", N)	#should try with conserve_qns=true per ed
	
	os = OpSum()
	for j in 1:(N-1)
		os += "Sz", j, "Sz", j+1
		os += 0.5, "S+", j, "S-", j+1
		os += 0.5, "S-", j, "S+", j+1
	end
	
	H = MPO(os, sites)
	
	psi0 = randomMPS(sites; linkdims=10)
	
	nsweeps = 10
	maxdim = [10,20,100,100,200]
	
	cutoff = [1e-6]
	
	energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
	
	#=commenting out some stuff hereafter
	@printf("Final energy = %.12f, n = %d\n", energy, N)
	
	open("heisenberg-one-half.txt","a") do io
   		@printf(io,"Final energy = %.12f, n = %d\n", energy, N)
	end
	=#
	energy
end

main(parse(Int64,ARGS[1]))
