using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(9999)

let
	N = 40
	J = 1		#see hamiltonian, should i use this
	h = 1		#see hamiltonian, should i use this
	
	sites = siteinds("S=1/2", N)
	
	os = OpSum()
	for j in 1:(N-1)
		os += "Sz", j, "Sz", j+1
		os += 0.5, "S+", j, "S-", j+1
		os += 0.5, "S-", j, "S+", j+1
	end
	
	H = MPO(os, sites)
	
	psi0 = randomMPS(sites; linkdims=30)	#linkdims is bond dimension
	
	nsweeps = 4 	#TODO why 5 initial? changed to four bc 5th ~ irrelevant
	maxdim = [10,20,100,200]	#maxdim is max bond dimension per sweep
	
	cutoff = [1E-11]		#truncation error goal of each sweep
	
	energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
	
	@printf("Final energy = %.12f\n", energy)	#TODO is this the precision you want
	
	open("heisenberg-one-half-n40.txt","a") do io	#a: write, create, append
   		@printf(io,"Final energy = %.12f\n", energy)
	end
end
