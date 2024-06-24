using ITensors, ITensorMPS
using Printf
using Random

Random.seed!(6969)

function dmrgLadder(numPairs, jl, jr, operator=-1)
	
	data = getOperator(numPairs, jl, jr, operator)
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

function getOperator(numPairs, jl, jr, operator)
	if operator == -1
		return hamiltonian(numPairs, jl, jr)
	end
	if typeof(operator) <: Array
		if true #(length(operator) <= 4)
			return siteOperator(numPairs, data)
		else
			return ArgumentError("operator unrecognised e#01")
		end
	end
	return ArgumentError("operator unrecognised e#02")
end

function hamiltonian(numPairs, Jl, Jr)
	N = Int(numPairs*2)

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

#TODO fix errors with tupling
function siteOperator(numPairs, opMatrix)
	operators1 = opMatrix[1]
	operators2 = opMatrix[2]	#TODO should i assert that operators1 and operators2 are the same length? probably :S
	siteSets1 = opMatrix[3]
	siteSets2 = opMatrix[4]
	
	sites = siteinds("S=1/2", numPairs*2)
	os = OpSum()
	numSiteSets = length(siteSets1)
	
	if numSiteSets < length(operators1)
		print("WARNING: one or more operator sets will be ignored as they have no corresponding site sets\n")
	end
	for i in 1:numSiteSets
		if isempty(operators2[1])
			operators2[i+1] = "I"
		end
		if isempty(operators1[i])
			operators1[i] = "I"
		end
		#get site op for this set
		for j in 1:length(siteSets1[i])	#TODO implement constant multiplication later
			print("hew")
			os += operators1[i], siteSets1[j], operators2[i], siteSets2[j]
		end	
	end
	
	op = MPO(os,sites)
	
	return [op, sites]
end

data = [("Sz",),("I",),(1,2),(1,2)]

#display(dmrgLadder(2,1,1,data)[1])
