using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(6969)

function dmrgLadder(length, jl, jr, operator=-1)
	
	data = getOperator(length, jl, jl, operator)
	op = data[1]
	sites = data[2]
	
	psi0 = randomMPS(sites; linkdims=10)
	
	nsweeps = 20#changed from 10
	maxdim = [200]#increased
	
	cutoff = [1e-12]#changed from 1e-6
	
	energy, psi = dmrg(op, psi0; nsweeps, maxdim, cutoff)
	
	#@printf("Final energy = %.12f, n = %d, Jl = %f, Jr = %f\n", energy, N, Jl, Jr)
	
	return [energy, psi]
end

function getOperator(length, jl, jr, operator)
	if operator == -1
		return hamiltonian(length, jl, jr)
	end
	if typeof(operator) <: Array
		if length(operator) >= 2
			return#TODO return something
		else
			return ArgumentError("operator unrecognised e#01")
		end
	end
	return ArgumentError("operator unrecognised e#02")
end

function hamiltonian(length, jl, jr)
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
	
	op = MPO(os, sites)
	
	return [op, sites]
end

function siteOperator(opMatrix)
	operators = opMatrix[1]
	opdSites = opMatrix[2]
	
	siteSets = length(opdSites)
	
	if siteSets < length(operators)
		print("WARNING: one or more operator sets will be ignored as they have no corresponding site sets\n")
	end
	for i in 1:siteSets
		#TODO apply op !	
	end
end

dmrgLadder(1,1,1)
